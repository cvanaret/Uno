// Detection of quadratic structure: every AMPL encoding of a quadratic must end up in the COO Hessian, with the same
// values and derivatives as the AD tapes (detect_quadratic_structure = false).
#include <map>
#include "test_utils.hpp"

using namespace nlwriter;

namespace {
   using Matrix = std::map<std::pair<int, int>, double>; // lower triangle of H (f = 1/2 x^T H x + ...)

   Matrix quadratic_form(const cppasl::QuadraticFormView& view) {
      Matrix matrix;
      for (std::size_t q = 0; q < view.number_nonzeros; ++q) matrix[{view.row_indices[q], view.column_indices[q]}] += view.values[q];
      return matrix;
   }

   bool same(const Matrix& a, const Matrix& b) {
      if (a.size() != b.size()) return false;
      for (const auto& [key, value]: a) {
         const auto other = b.find(key);
         if (other == b.end() || !testing::close(value, other->second, 1e-14)) return false;
      }
      return true;
   }

   struct Case {
      const char* name;
      Expression expression;
      Matrix hessian;
      double constant_at_zero;
   };

   std::vector<Case> encodings() {
      return {
         {"product of variables", num(0.5) * op_list(54, {num(2.) * var(0) * var(0), num(3.) * var(0) * var(1), num(3.) * var(1) * var(0)}),
            {{{0, 0}, 2.}, {{1, 0}, 3.}}, 0.},
         {"power 2", pow(var(0), num(2.)), {{{0, 0}, 2.}}, 0.},
         {"square of an affine form", pow(var(0) + num(2.) * var(1) - num(1.), num(2.)),
            {{{0, 0}, 2.}, {{1, 0}, 4.}, {{1, 1}, 8.}}, 1.},
         {"product of affine forms", (var(0) + var(1)) * (var(1) - num(3.)), {{{1, 0}, 1.}, {{1, 1}, 2.}}, 0.},
         {"division by a constant", var(0) * var(1) / num(4.), {{{1, 0}, 0.25}}, 0.},
         {"negation", neg(var(0) * var(0)), {{{0, 0}, -2.}}, 0.},
         {"merged duplicates", num(3.) * (var(0) * var(1)) - var(1) * var(0), {{{1, 0}, 2.}}, 0.},
         {"cancellation", var(0) * var(1) - var(1) * var(0) + var(2), {}, 0.},
         {"defined variable", var(3) * var(3) + var(4), {{{0, 0}, 2.}, {{1, 0}, 2.}, {{1, 1}, 2.}}, 1.},
         {"power 1 and constant exponents", pow(var(2), num(1.)) * pow(num(2.), num(2.)) * var(1), {{{2, 1}, 4.}}, 0.},
         {"nested sums and scalings", num(2.) * op_list(54, {var(0), var(1) * var(2), neg(var(2) * num(-1.) * var(2))}),
            {{{2, 1}, 2.}, {{2, 2}, 4.}}, 0.},
      };
   }

   NlProblem problem_with(const Expression& expression) {
      NlProblem p;
      p.number_variables = 3;
      p.defined_variables.push_back({{{0, 1.}, {1, 1.}}, {}});   // v3 = x0 + x1
      p.defined_variables.push_back({{}, num(1.)});               // v4 = 1 (constant)
      p.objectives.push_back({expression, {}});
      p.objective_senses = {0};
      p.constraints.push_back({expression, {{2, 1.}}}); // the same function as a constraint
      p.constraint_lower = {-INFINITY};
      p.constraint_upper = {1.};
      return p;
   }
} // namespace

TEST_CASE(quadratic_encodings_are_detected) {
   std::mt19937_64 random(3);
   for (const Case& c: encodings()) {
      const NlProblem problem = problem_with(c.expression);
      const cppasl::NlModel model = testing::write_and_read(problem, Format::Text, "quadratic");
      cppasl::ReaderOptions tape_options;
      tape_options.detect_quadratic_structure = false;
      const cppasl::NlModel tapes = testing::write_and_read(problem, Format::Binary, "quadratic", tape_options);
      if (!CHECK(same(quadratic_form(model.objective_quadratic_form()), c.hessian))) std::fprintf(stderr, "    case: %s\n", c.name);
      CHECK(same(quadratic_form(model.constraint_quadratic_form(0)), c.hessian));
      CHECK(model.objective_structure() == (c.hessian.empty() ? cppasl::FunctionStructure::Linear : cppasl::FunctionStructure::Quadratic));
      CHECK(tapes.objective_quadratic_form().number_nonzeros == 0);
      cppasl::EvaluationWorkspace w1(model), w2(tapes);
      const std::vector<double> zero(3, 0.);
      CHECK_CLOSE(model.evaluate_objective(w1, zero.data()), c.constant_at_zero, 1e-14);
      for (int trial = 0; trial < 3; ++trial) {
         const std::vector<double> x = testing::random_vector(random, 3, -2., 2.);
         CHECK_CLOSE(model.evaluate_objective(w1, x.data()), tapes.evaluate_objective(w2, x.data()), 1e-13);
         std::vector<double> g1(3), g2(3), c1(1), c2(1);
         model.evaluate_objective_gradient(w1, x.data(), g1.data());
         tapes.evaluate_objective_gradient(w2, x.data(), g2.data());
         for (int j = 0; j < 3; ++j) CHECK_CLOSE(g1[static_cast<std::size_t>(j)], g2[static_cast<std::size_t>(j)], 1e-13);
         model.evaluate_constraints(w1, x.data(), c1.data());
         tapes.evaluate_constraints(w2, x.data(), c2.data());
         CHECK_CLOSE(c1[0], c2[0], 1e-13);
         const std::vector<double> y{-1.5};
         std::vector<double> h1(model.number_hessian_nonzeros()), h2(tapes.number_hessian_nonzeros());
         model.evaluate_lagrangian_hessian(w1, x.data(), 0.5, y.data(), h1.data());
         tapes.evaluate_lagrangian_hessian(w2, x.data(), 0.5, y.data(), h2.data());
         const std::vector<double> d1 = testing::dense_hessian(model, h1), d2 = testing::dense_hessian(tapes, h2);
         for (std::size_t k = 0; k < 9; ++k) CHECK_CLOSE(d1[k], d2[k], 1e-13);
         CHECK_DERIVATIVES(model, x, 0.5, y, 1e-6);
         CHECK_DERIVATIVES(tapes, x, 0.5, y, 1e-6);
      }
   }
}

TEST_CASE(quadratic_form_layout) {
   // diagonal first, then strictly lower sorted by (column, row)
   const NlProblem problem = problem_with(op_list(54, {var(2) * var(0), var(1) * var(1), var(1) * var(0), var(2) * var(2), var(2) * var(1)}));
   const cppasl::NlModel model = testing::write_and_read(problem, Format::Text, "layout");
   const cppasl::QuadraticFormView view = model.objective_quadratic_form();
   CHECK(view.number_nonzeros == 5 && view.number_diagonal_nonzeros == 2);
   const std::vector<std::pair<int, int>> expected{{1, 1}, {2, 2}, {1, 0}, {2, 0}, {2, 1}};
   for (std::size_t q = 0; q < 5; ++q) CHECK(std::make_pair(view.row_indices[q], view.column_indices[q]) == expected[q]);
   CHECK(view.values[0] == 2. && view.values[2] == 1.);
}

TEST_CASE(mixed_quadratic_and_nonlinear_terms) {
   // f = x0^2 + x0 x1 + exp(x1) + 3 x2 : quadratic part in COO, exp in one element, linear part merged
   const NlProblem problem = problem_with(op_list(54, {var(0) * var(0), var(0) * var(1), op(44, {var(1)}), num(3.) * var(2)}));
   const cppasl::NlModel model = testing::write_and_read(problem, Format::Text, "mixed");
   CHECK(model.objective_structure() == cppasl::FunctionStructure::Nonlinear);
   CHECK(model.objective_quadratic_form().number_nonzeros == 2);
   std::mt19937_64 random(4);
   const std::vector<double> x = testing::random_vector(random, 3);
   cppasl::EvaluationWorkspace workspace(model);
   CHECK_CLOSE(model.evaluate_objective(workspace, x.data()), x[0] * x[0] + x[0] * x[1] + std::exp(x[1]) + 3. * x[2], 1e-14);
   CHECK_DERIVATIVES(model, x, 1., std::vector<double>{2.}, 1e-6);
}
