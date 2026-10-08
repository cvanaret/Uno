// Large random instances: the reader and the evaluations scale linearly and agree with direct computations.
#include <chrono>
#include "test_utils.hpp"

using namespace nlwriter;

namespace {
   double seconds_since(std::chrono::steady_clock::time_point start) {
      return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
   }

   struct Triplet {
      int row, column;
      double value;
   };

   /// 1/2 x^T Q x encoded as AMPL does: o2 n0.5 o54 {o2 o2 n v v} (both triangles)
   Expression encode(const std::vector<Triplet>& q) {
      std::vector<Expression> terms;
      terms.reserve(2 * q.size());
      for (const Triplet& t: q) {
         terms.push_back(num(t.value) * var(t.row) * var(t.column));
         if (t.row != t.column) terms.push_back(num(t.value) * var(t.column) * var(t.row));
      }
      return num(0.5) * op_list(54, terms);
   }

   /// random symmetric Q (lower triplets), diagonally dominant
   std::vector<Triplet> random_matrix(std::mt19937_64& random, int n, int off_diagonal_per_row) {
      std::uniform_int_distribution<int> column(0, n - 1);
      std::uniform_real_distribution<double> value(-1., 1.);
      std::vector<Triplet> q;
      std::vector<double> absolute_sums(static_cast<std::size_t>(n), 0.);
      for (int i = 0; i < n; ++i) {
         for (int k = 0; k < off_diagonal_per_row; ++k) {
            const int j = column(random);
            if (j >= i) continue;
            q.push_back({i, j, value(random)});
            absolute_sums[static_cast<std::size_t>(i)] += std::fabs(q.back().value);
            absolute_sums[static_cast<std::size_t>(j)] += std::fabs(q.back().value);
         }
      }
      for (int i = 0; i < n; ++i) q.push_back({i, i, absolute_sums[static_cast<std::size_t>(i)] + 0.5});
      return q;
   }

   std::vector<double> multiply(const std::vector<Triplet>& q, const std::vector<double>& x) {
      std::vector<double> y(x.size(), 0.);
      for (const Triplet& t: q) {
         y[static_cast<std::size_t>(t.row)] += t.value * x[static_cast<std::size_t>(t.column)];
         if (t.row != t.column) y[static_cast<std::size_t>(t.column)] += t.value * x[static_cast<std::size_t>(t.row)];
      }
      return y;
   }
} // namespace

TEST_CASE(large_random_qp_and_qcqp) {
   for (int n: {1000, 20000, 100000}) {
      std::mt19937_64 random(static_cast<std::uint64_t>(n));
      const std::vector<Triplet> Q = random_matrix(random, n, 4);
      const int number_quadratic_constraints = 5, number_linear_constraints = n / 10;
      std::vector<std::vector<Triplet>> P;
      NlProblem p;
      p.number_variables = n;
      p.objectives.push_back({encode(Q), {{0, 1.}, {n - 1, -2.}}});
      std::uniform_int_distribution<int> column(0, n - 1);
      for (int i = 0; i < number_quadratic_constraints; ++i) {
         P.push_back(random_matrix(random, n, 2));
         p.constraints.push_back({encode(P.back()), {{column(random), 1.}}});
         p.constraint_lower.push_back(-INFINITY);
         p.constraint_upper.push_back(1.);
      }
      for (int i = 0; i < number_linear_constraints; ++i) {
         p.constraints.push_back({{}, {{column(random), 1.}, {column(random), -1.}, {column(random), 0.5}}});
         p.constraint_lower.push_back(-1.);
         p.constraint_upper.push_back(1.);
      }
      const std::string contents = p.write(Format::Text);
      const auto start = std::chrono::steady_clock::now();
      const cppasl::NlModel model = cppasl::NlModel::read_from_memory(contents.data(), contents.size());
      const double read_time = seconds_since(start);
      const auto classify_start = std::chrono::steady_clock::now();
      const std::string type = model.classify().to_string();
      const double classify_time = seconds_since(classify_start);
      std::printf("    n = %6d: %8.1f MB read in %.3f s, classified \"%s\" in %.3f s\n", n,
         static_cast<double>(contents.size()) / 1e6, read_time, type.c_str(), classify_time);
      CHECK(type == "convex QCQP");

      cppasl::EvaluationWorkspace workspace(model);
      const std::vector<double> x = testing::random_vector(random, static_cast<std::size_t>(n));
      const std::vector<double> Qx = multiply(Q, x);
      double expected = x[0] - 2. * x[static_cast<std::size_t>(n - 1)];
      for (int j = 0; j < n; ++j) expected += 0.5 * x[static_cast<std::size_t>(j)] * Qx[static_cast<std::size_t>(j)];
      CHECK_CLOSE(model.evaluate_objective(workspace, x.data()), expected, 1e-12);
      std::vector<double> gradient(static_cast<std::size_t>(n));
      model.evaluate_objective_gradient(workspace, x.data(), gradient.data());
      double gradient_error = 0.;
      for (int j = 0; j < n; ++j) {
         const double reference = Qx[static_cast<std::size_t>(j)] + (j == 0 ? 1. : 0.) + (j == n - 1 ? -2. : 0.);
         gradient_error = std::max(gradient_error, std::fabs(gradient[static_cast<std::size_t>(j)] - reference));
      }
      CHECK(gradient_error < 1e-10);
      // constraint values, and HVP of the Lagrangian against the direct sparse products
      std::vector<double> constraints(static_cast<std::size_t>(model.number_constraints()));
      model.evaluate_constraints(workspace, x.data(), constraints.data());
      const std::vector<double> y = testing::random_vector(random, constraints.size());
      const std::vector<double> v = testing::random_vector(random, static_cast<std::size_t>(n));
      std::vector<double> reference = multiply(Q, v);
      for (int i = 0; i < number_quadratic_constraints; ++i) {
         const std::vector<double> Px = multiply(P[static_cast<std::size_t>(i)], x), Pv = multiply(P[static_cast<std::size_t>(i)], v);
         double quadratic = 0.;
         for (int j = 0; j < n; ++j) {
            quadratic += 0.5 * x[static_cast<std::size_t>(j)] * Px[static_cast<std::size_t>(j)];
            reference[static_cast<std::size_t>(j)] += y[static_cast<std::size_t>(i)] * Pv[static_cast<std::size_t>(j)];
         }
         const double linear = constraints[static_cast<std::size_t>(i)] - quadratic;
         CHECK(std::fabs(linear) <= 1. + 1e-9); // one variable with coefficient 1
      }
      std::vector<double> hvp(static_cast<std::size_t>(n)), wv(static_cast<std::size_t>(n)), h(model.number_hessian_nonzeros());
      model.evaluate_lagrangian_hessian_vector_product(workspace, x.data(), 1., y.data(), v.data(), hvp.data());
      model.evaluate_lagrangian_hessian(workspace, x.data(), 1., y.data(), h.data());
      model.multiply_hessian(h.data(), v.data(), wv.data());
      double hvp_error = 0.;
      for (std::size_t j = 0; j < hvp.size(); ++j) {
         hvp_error = std::max({hvp_error, std::fabs(hvp[j] - reference[j]), std::fabs(wv[j] - reference[j])});
      }
      CHECK(hvp_error < 1e-10);
      if (n <= 1000) CHECK_DERIVATIVES(model, x, 1., y, 1e-5);
   }
}

TEST_CASE(large_partially_separable_nlp) {
   // f = sum_i exp(x_i - x_{i+1}) + sum_i x_i^4 ; c_i = sin(x_i) * x_{i+1}^2 for i < n/2
   const int n = 100000;
   NlProblem p;
   p.number_variables = n;
   std::vector<Expression> terms;
   for (int i = 0; i + 1 < n; ++i) terms.push_back(op(44, {var(i) - var(i + 1)}));
   for (int i = 0; i < n; ++i) terms.push_back(pow(var(i), num(4.)));
   p.objectives.push_back({op_list(54, terms), {}});
   for (int i = 0; i < n / 2; ++i) {
      p.constraints.push_back({op(41, {var(i)}) * pow(var(i + 1), num(2.)), {}});
      p.constraint_lower.push_back(-1.);
      p.constraint_upper.push_back(1.);
   }
   const std::string contents = p.write(Format::Binary);
   const auto start = std::chrono::steady_clock::now();
   const cppasl::NlModel model = cppasl::NlModel::read_from_memory(contents.data(), contents.size());
   std::printf("    %.1f MB read in %.3f s, %zu Hessian nonzeros\n", static_cast<double>(contents.size()) / 1e6,
      seconds_since(start), model.number_hessian_nonzeros());
   CHECK(model.number_hessian_nonzeros() == static_cast<std::size_t>(2 * n - 1));
   CHECK(model.classify().to_string() == "NLP");
   std::mt19937_64 random(21);
   const std::vector<double> x = testing::random_vector(random, n);
   cppasl::EvaluationWorkspace workspace(model);
   double expected = 0.;
   for (int i = 0; i + 1 < n; ++i) expected += std::exp(x[static_cast<std::size_t>(i)] - x[static_cast<std::size_t>(i + 1)]);
   for (int i = 0; i < n; ++i) expected += std::pow(x[static_cast<std::size_t>(i)], 4);
   CHECK_CLOSE(model.evaluate_objective(workspace, x.data()), expected, 1e-12);
   const std::vector<double> y = testing::random_vector(random, static_cast<std::size_t>(n / 2));
   std::vector<double> h(model.number_hessian_nonzeros()), hvp(static_cast<std::size_t>(n)), wv(static_cast<std::size_t>(n));
   const std::vector<double> v = testing::random_vector(random, static_cast<std::size_t>(n));
   const auto hessian_start = std::chrono::steady_clock::now();
   model.evaluate_lagrangian_hessian(workspace, x.data(), 1., y.data(), h.data());
   const double hessian_time = seconds_since(hessian_start);
   model.evaluate_lagrangian_hessian_vector_product(workspace, x.data(), 1., y.data(), v.data(), hvp.data());
   model.multiply_hessian(h.data(), v.data(), wv.data());
   double error = 0.;
   for (std::size_t j = 0; j < hvp.size(); ++j) error = std::max(error, std::fabs(hvp[j] - wv[j]) / std::max(1., std::fabs(wv[j])));
   CHECK(error < 1e-12);
   std::printf("    Lagrangian Hessian in %.3f s\n", hessian_time);
}
