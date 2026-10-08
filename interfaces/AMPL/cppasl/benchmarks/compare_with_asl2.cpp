// Cross-validation and performance comparison of cppasl against the original ASL2 (solvers2, pfgh_read).
// usage: compare_with_asl2 file.nl [repetitions] [--quiet]
// Values are compared entry by entry through (row, column) maps; timings are the best of the repetitions,
// alternating between two points so that no implementation can reuse cached values.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <random>
#include <string>
#include <vector>
#include "asl2_session.h"
#include "cppasl/nl_model.hpp"

namespace {
   using Clock = std::chrono::steady_clock;
   double seconds(Clock::time_point start) { return std::chrono::duration<double>(Clock::now() - start).count(); }

   struct Comparison {
      double worst{0.};
      std::string where;
      void add(double a, double b, const std::string& what) {
         const double error = std::fabs(a - b) / std::max(1., std::max(std::fabs(a), std::fabs(b)));
         if (!(error <= this->worst)) { // also catches NaN
            this->worst = std::isnan(error) ? INFINITY : error;
            this->where = what;
         }
      }
   };

   using SparseMap = std::map<std::pair<int, int>, double>;

   /// fastest of `repetitions` calls to body(k) (k alternates between the two points), prepare(k) untimed before
   double best_time(int repetitions, const std::function<void(int)>& body, const std::function<void(int)>& prepare = {}) {
      double best = INFINITY;
      for (int r = 0; r < repetitions; ++r) {
         if (prepare) prepare(r % 2);
         const auto start = Clock::now();
         body(r % 2);
         best = std::min(best, seconds(start));
      }
      return best;
   }
} // namespace

int main(int argc, char** argv) {
   if (argc < 2) {
      std::fprintf(stderr, "usage: %s file.nl [repetitions] [--quiet]\n", argv[0]);
      return 1;
   }
   const std::string file = argv[1];
   const int repetitions = (argc > 2) ? std::max(1, std::atoi(argv[2])) : 5;
   const bool quiet = (argc > 3) && std::strcmp(argv[3], "--quiet") == 0;

   // ---------------------------------------------------------------------------------------------- reading
   auto start = Clock::now();
   Asl2Session* asl = asl2_read(file.c_str());
   const double asl_read = seconds(start);
   if (asl == nullptr) { std::fprintf(stderr, "ASL2 cannot read %s\n", file.c_str()); return 1; }
   start = Clock::now();
   const std::size_t asl_hessian_nonzeros = asl2_hessian_setup(asl);
   const double asl_hessian_setup = seconds(start);

   start = Clock::now();
   const cppasl::NlModel model = cppasl::NlModel::read_file(file);
   const double cppasl_read = seconds(start);
   start = Clock::now();
   const cppasl::ProblemClassification classification = model.classify();
   const double cppasl_classify = seconds(start);
   cppasl::EvaluationWorkspace workspace(model);

   const int n = model.number_variables(), m = model.number_constraints();
   const auto un = static_cast<std::size_t>(n), um = static_cast<std::size_t>(m);
   if (asl2_number_variables(asl) != n || asl2_number_constraints(asl) != m) {
      std::fprintf(stderr, "dimension mismatch\n");
      return 1;
   }

   // two evaluation points around the initial point (kept inside (-1, 1) for bounded domains), random multipliers
   std::mt19937_64 random(2024);
   std::uniform_real_distribution<double> perturbation(-0.5, 0.5), multiplier(-1., 1.);
   std::vector<double> x0(un);
   asl2_initial_point(asl, x0.data());
   std::vector<std::vector<double>> points(2, x0), multipliers(2, std::vector<double>(um));
   for (auto& point: points) for (double& value: point) value = (std::fabs(value) < 1e-12) ? 0.8 * perturbation(random) : value * (1. + 0.1 * perturbation(random));
   for (auto& y: multipliers) for (double& value: y) value = multiplier(random);
   const double sigma[2] = {1., 0.75};
   std::vector<double> direction(un);
   for (double& value: direction) value = perturbation(random);

   // ---------------------------------------------------------------------------------------------- correctness
   Comparison objective, gradient, constraints, jacobian, hessian, hvp;
   const std::size_t asl_jacobian_nonzeros = asl2_number_jacobian_nonzeros(asl);
   std::vector<int> asl_jacobian_rows(asl_jacobian_nonzeros), asl_jacobian_columns(asl_jacobian_nonzeros);
   asl2_jacobian_structure(asl, asl_jacobian_rows.data(), asl_jacobian_columns.data());
   std::vector<int> asl_hessian_rows(asl_hessian_nonzeros), asl_hessian_columns(asl_hessian_nonzeros);
   asl2_hessian_structure(asl, asl_hessian_rows.data(), asl_hessian_columns.data());

   std::vector<double> ga(un), gc(un), ca(um), cc(um), ja(asl_jacobian_nonzeros), jc(model.number_jacobian_nonzeros());
   std::vector<double> ha(asl_hessian_nonzeros), hc(model.number_hessian_nonzeros()), va(un), vc(un);
   auto compare_sparse = [](Comparison& comparison, const SparseMap& reference, const SparseMap& ours, const char* name) {
      SparseMap all(reference);
      for (const auto& entry: ours) all.emplace(entry.first, 0.);
      for (const auto& [key, value]: all) {
         const auto other = ours.find(key);
         comparison.add(value, other == ours.end() ? 0. : other->second,
            std::string(name) + "(" + std::to_string(key.first) + "," + std::to_string(key.second) + ")");
      }
   };
   for (int k = 0; k < 2; ++k) {
      const double* x = points[static_cast<std::size_t>(k)].data();
      const double* y = multipliers[static_cast<std::size_t>(k)].data();
      if (model.number_objectives() > 0) {
         objective.add(asl2_objective(asl, x), model.evaluate_objective(workspace, x), "f");
         asl2_objective_gradient(asl, x, ga.data());
         model.evaluate_objective_gradient(workspace, x, gc.data());
         for (std::size_t j = 0; j < un; ++j) gradient.add(ga[j], gc[j], "g[" + std::to_string(j) + "]");
      }
      if (m > 0) {
         asl2_constraints(asl, x, ca.data());
         model.evaluate_constraints(workspace, x, cc.data());
         for (std::size_t i = 0; i < um; ++i) constraints.add(ca[i], cc[i], "c[" + std::to_string(i) + "]");
         asl2_jacobian(asl, x, ja.data());
         model.evaluate_jacobian(workspace, x, jc.data());
         SparseMap asl_map, our_map;
         for (std::size_t q = 0; q < ja.size(); ++q) asl_map[{asl_jacobian_rows[q], asl_jacobian_columns[q]}] += ja[q];
         for (std::size_t q = 0; q < jc.size(); ++q) our_map[{model.jacobian_row_indices()[q], model.jacobian_column_indices()[q]}] += jc[q];
         compare_sparse(jacobian, asl_map, our_map, "J");
      }
      // Hessian of the Lagrangian (ASL: upper triangle; cppasl: lower triangle)
      if (model.number_objectives() > 0) (void)asl2_objective(asl, x);
      if (m > 0) asl2_constraints(asl, x, ca.data());
      asl2_hessian(asl, x, sigma[k], y, ha.data());
      model.evaluate_lagrangian_hessian(workspace, x, sigma[k], y, hc.data());
      SparseMap asl_map, our_map;
      for (std::size_t q = 0; q < ha.size(); ++q) asl_map[{asl_hessian_columns[q], asl_hessian_rows[q]}] += ha[q];
      for (std::size_t q = 0; q < hc.size(); ++q) our_map[{model.hessian_row_indices()[q], model.hessian_column_indices()[q]}] += hc[q];
      compare_sparse(hessian, asl_map, our_map, "H");
      if (model.number_objectives() > 0) (void)asl2_objective(asl, x);
      if (m > 0) asl2_constraints(asl, x, ca.data());
      asl2_hessian_vector_product(asl, x, sigma[k], y, direction.data(), va.data());
      model.evaluate_lagrangian_hessian_vector_product(workspace, x, sigma[k], y, direction.data(), vc.data());
      for (std::size_t j = 0; j < un; ++j) hvp.add(va[j], vc[j], "Hv[" + std::to_string(j) + "]");
   }

   // ---------------------------------------------------------------------------------------------- performance
   struct Row {
      const char* name;
      double asl, ours;
      const Comparison* comparison;
   };
   std::vector<Row> rows;
   rows.push_back({"read", asl_read, cppasl_read, nullptr});
   rows.push_back({"Hessian structure", asl_hessian_setup, 0., nullptr});
   auto point = [&](int k) { return points[static_cast<std::size_t>(k)].data(); };
   auto y_of = [&](int k) { return multipliers[static_cast<std::size_t>(k)].data(); };
   auto asl_prepare = [&](int k) { // ASL differentiates at the point of the last function evaluation
      if (model.number_objectives() > 0) (void)asl2_objective(asl, point(k));
      if (m > 0) asl2_constraints(asl, point(k), ca.data());
   };
   auto our_prepare = [&](int k) {
      if (model.number_objectives() > 0) (void)model.evaluate_objective(workspace, point(k));
      if (m > 0) model.evaluate_constraints(workspace, point(k), cc.data());
   };
   const int r = 2 * repetitions;
   if (model.number_objectives() > 0) {
      rows.push_back({"objective", best_time(r, [&](int k) { (void)asl2_objective(asl, point(k)); }),
         best_time(r, [&](int k) { (void)model.evaluate_objective(workspace, point(k)); }), &objective});
      rows.push_back({"objective gradient", best_time(r, [&](int k) { asl2_objective_gradient(asl, point(k), ga.data()); }, asl_prepare),
         best_time(r, [&](int k) { model.evaluate_objective_gradient(workspace, point(k), gc.data()); }, our_prepare), &gradient});
   }
   if (m > 0) {
      rows.push_back({"constraints", best_time(r, [&](int k) { asl2_constraints(asl, point(k), ca.data()); }),
         best_time(r, [&](int k) { model.evaluate_constraints(workspace, point(k), cc.data()); }), &constraints});
      rows.push_back({"Jacobian", best_time(r, [&](int k) { asl2_jacobian(asl, point(k), ja.data()); }, asl_prepare),
         best_time(r, [&](int k) { model.evaluate_jacobian(workspace, point(k), jc.data()); }, our_prepare), &jacobian});
   }
   rows.push_back({"Lagrangian Hessian", best_time(r, [&](int k) { asl2_hessian(asl, point(k), sigma[k], y_of(k), ha.data()); }, asl_prepare),
      best_time(r, [&](int k) { model.evaluate_lagrangian_hessian(workspace, point(k), sigma[k], y_of(k), hc.data()); }, our_prepare), &hessian});
   rows.push_back({"Hessian-vector product", best_time(r, [&](int k) {
         asl2_hessian_vector_product(asl, point(k), sigma[k], y_of(k), direction.data(), va.data()); }, asl_prepare),
      best_time(r, [&](int k) {
         model.evaluate_lagrangian_hessian_vector_product(workspace, point(k), sigma[k], y_of(k), direction.data(), vc.data()); }, our_prepare), &hvp});

   std::printf("%s: n = %d, m = %d, Jacobian nnz %zu (ASL %zu), Hessian nnz %zu (ASL %zu), %s (classified in %.3f s)\n",
      file.c_str(), n, m, model.number_jacobian_nonzeros(), asl_jacobian_nonzeros, model.number_hessian_nonzeros(),
      asl_hessian_nonzeros, classification.to_string().c_str(), cppasl_classify);
   double worst = 0.;
   if (!quiet) std::printf("  %-24s %12s %12s %9s %11s\n", "operation", "ASL2 [s]", "cppasl [s]", "speedup", "max rel err");
   for (const Row& row: rows) {
      const double error = row.comparison ? row.comparison->worst : 0.;
      worst = std::max(worst, error);
      if (quiet) continue;
      if (row.ours > 0.) std::printf("  %-24s %12.6f %12.6f %8.1fx %11.2e\n", row.name, row.asl, row.ours, row.asl / row.ours, error);
      else std::printf("  %-24s %12.6f %12s %9s %11s\n", row.name, row.asl, "(in read)", "", "");
   }
   for (const Row& row: rows) {
      if (row.comparison && row.comparison->worst > 1e-9) std::printf("  largest %s discrepancy at %s\n", row.name, row.comparison->where.c_str());
   }
   std::printf("  => %s (max relative difference %.2e)\n", worst <= 1e-9 ? "AGREE" : "DISAGREE", worst);
   asl2_close(asl);
   return worst <= 1e-9 ? 0 : 2;
}
