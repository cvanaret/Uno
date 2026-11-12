#include <chrono>
#include <cstring>
#include <exception>
#include <fstream>
#include "test_utils.hpp"

namespace testing {
   std::vector<TestCase>& registry() {
      static std::vector<TestCase> tests;
      return tests;
   }
   int failed_checks = 0;
   int total_checks = 0;

   cppasl::NlModel write_and_read(const nlwriter::NlProblem& problem, nlwriter::Format format, const std::string& name,
         const cppasl::ReaderOptions& options) {
      const std::string contents = problem.write(format);
      if (const char* directory = std::getenv("CPPASL_TEST_OUTPUT")) {
         nlwriter::write_file(std::string(directory) + "/" + name + (format == nlwriter::Format::Text ? "" : "_binary") + ".nl", contents);
      }
      return cppasl::NlModel::read_from_memory(contents.data(), contents.size(), options);
   }

   std::vector<double> dense_hessian(const cppasl::NlModel& model, const std::vector<double>& values) {
      const auto n = static_cast<std::size_t>(model.number_variables());
      std::vector<double> dense(n * n, 0.);
      for (std::size_t k = 0; k < values.size(); ++k) {
         const auto i = static_cast<std::size_t>(model.hessian_row_indices()[k]);
         const auto j = static_cast<std::size_t>(model.hessian_column_indices()[k]);
         dense[i * n + j] += values[k];
         if (i != j) dense[j * n + i] += values[k];
      }
      return dense;
   }

   double check_derivatives(const cppasl::NlModel& model, const std::vector<double>& x0, double sigma,
         const std::vector<double>& y, double tolerance, const char* file, int line) {
      const auto n = static_cast<std::size_t>(model.number_variables());
      const auto m = static_cast<std::size_t>(model.number_constraints());
      cppasl::EvaluationWorkspace workspace(model);
      std::vector<double> x(x0), gradient(n), constraints(m), jacobian(model.number_jacobian_nonzeros());
      double worst = 0.;
      auto compare = [&](double computed, double reference, const char* what, std::size_t i, std::size_t j) {
         const double error = std::fabs(computed - reference) / std::max(1., std::fabs(reference));
         worst = std::max(worst, error);
         if (error > tolerance) {
            report(false, what, file, line, (std::ostringstream() << "[" << i << "," << j << "] " << computed << " vs " << reference).str());
            return false;
         }
         return true;
      };
      auto lagrangian_gradient = [&](const std::vector<double>& point) {
         std::vector<double> g(n), jv(model.number_jacobian_nonzeros()), jty(n);
         if (model.number_objectives() > 0) model.evaluate_objective_gradient(workspace, point.data(), g.data());
         for (double& value: g) value *= sigma;
         if (m > 0) {
            model.evaluate_jacobian(workspace, point.data(), jv.data());
            model.multiply_jacobian_transpose(jv.data(), y.data(), jty.data());
            for (std::size_t i = 0; i < n; ++i) g[i] += jty[i];
         }
         return g;
      };
      if (model.number_objectives() > 0) model.evaluate_objective_gradient(workspace, x.data(), gradient.data());
      if (m > 0) model.evaluate_jacobian(workspace, x.data(), jacobian.data());
      std::vector<double> dense_jacobian(m * n, 0.);
      for (std::size_t k = 0; k < jacobian.size(); ++k) {
         dense_jacobian[static_cast<std::size_t>(model.jacobian_row_indices()[k]) * n + static_cast<std::size_t>(model.jacobian_column_indices()[k])] += jacobian[k];
      }
      std::vector<double> hessian_values(model.number_hessian_nonzeros());
      model.evaluate_lagrangian_hessian(workspace, x.data(), sigma, y.data(), hessian_values.data());
      const std::vector<double> hessian = dense_hessian(model, hessian_values);
      const std::vector<double> reference_gradient = lagrangian_gradient(x);

      std::vector<double> plus(m), minus(m);
      for (std::size_t j = 0; j < n; ++j) {
         const double h = 1e-6 * std::max(1., std::fabs(x0[j]));
         x[j] = x0[j] + h;
         const double f_plus = model.number_objectives() > 0 ? model.evaluate_objective(workspace, x.data()) : 0.;
         if (m > 0) model.evaluate_constraints(workspace, x.data(), plus.data());
         const std::vector<double> g_plus = lagrangian_gradient(x);
         x[j] = x0[j] - h;
         const double f_minus = model.number_objectives() > 0 ? model.evaluate_objective(workspace, x.data()) : 0.;
         if (m > 0) model.evaluate_constraints(workspace, x.data(), minus.data());
         const std::vector<double> g_minus = lagrangian_gradient(x);
         x[j] = x0[j];
         if (model.number_objectives() > 0) compare(gradient[j], (f_plus - f_minus) / (2. * h), "objective gradient", j, 0);
         for (std::size_t i = 0; i < m; ++i) compare(dense_jacobian[i * n + j], (plus[i] - minus[i]) / (2. * h), "Jacobian", i, j);
         for (std::size_t i = 0; i < n; ++i) compare(hessian[i * n + j], (g_plus[i] - g_minus[i]) / (2. * h), "Hessian", i, j);
      }
      (void)reference_gradient;
      // HVP (matrix-free) versus the assembled Hessian times v
      std::mt19937_64 random(7);
      const std::vector<double> v = random_vector(random, n);
      std::vector<double> hvp(n), wv(n);
      model.evaluate_lagrangian_hessian_vector_product(workspace, x.data(), sigma, y.data(), v.data(), hvp.data());
      model.multiply_hessian(hessian_values.data(), v.data(), wv.data());
      for (std::size_t i = 0; i < n; ++i) compare(hvp[i], wv[i], "Hessian-vector product", i, 0);
      // J v versus dense
      if (m > 0) {
         std::vector<double> jv(m);
         model.multiply_jacobian(jacobian.data(), v.data(), jv.data());
         for (std::size_t i = 0; i < m; ++i) {
            double reference = 0.;
            for (std::size_t j = 0; j < n; ++j) reference += dense_jacobian[i * n + j] * v[j];
            compare(jv[i], reference, "J v", i, 0);
         }
      }
      report(true, "derivatives", file, line);
      return worst;
   }
} // namespace testing

int main(int argc, char** argv) {
   const char* filter = (argc > 1) ? argv[1] : nullptr;
   int failed_tests = 0, run = 0;
   for (const testing::TestCase& test: testing::registry()) {
      if (filter != nullptr && std::strstr(test.name, filter) == nullptr) continue;
      const int failed_before = testing::failed_checks;
      const auto start = std::chrono::steady_clock::now();
      try {
         test.body();
      }
      catch (const std::exception& exception) {
         ++testing::failed_checks;
         std::fprintf(stderr, "  EXCEPTION: %s\n", exception.what());
      }
      const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
      const bool passed = testing::failed_checks == failed_before;
      failed_tests += !passed;
      ++run;
      std::printf("[%s] %-45s %8.3f s\n", passed ? "  OK  " : "FAILED", test.name, seconds);
   }
   std::printf("\n%d/%d tests passed, %d checks, %d failed\n", run - failed_tests, run, testing::total_checks, testing::failed_checks);
   return failed_tests == 0 ? 0 : 1;
}
