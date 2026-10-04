// First and second derivatives of every supported operation, checked by finite differences, plus tape reuse and
// defined variables shared by several functions.
#include "test_utils.hpp"

using namespace nlwriter;

namespace {
   struct UnaryCase {
      int opcode;
      const char* name;
      double low, high; // domain of the argument
   };
   const UnaryCase unary_cases[] = {
      {37, "tanh", -2., 2.}, {38, "tan", -1., 1.}, {39, "sqrt", 0.5, 3.}, {40, "sinh", -2., 2.}, {41, "sin", -3., 3.},
      {42, "log10", 0.5, 3.}, {43, "log", 0.5, 3.}, {44, "exp", -2., 2.}, {45, "cosh", -2., 2.}, {46, "cos", -3., 3.},
      {47, "atanh", -0.8, 0.8}, {49, "atan", -2., 2.}, {50, "asinh", -2., 2.}, {51, "asin", -0.8, 0.8},
      {52, "acosh", 1.5, 3.}, {53, "acos", -0.8, 0.8}, {15, "abs", 0.3, 2.}, {79, "logistic", -2., 2.}, {16, "neg", -2., 2.}};
} // namespace

TEST_CASE(unary_operations) {
   std::mt19937_64 random(11);
   for (const UnaryCase& c: unary_cases) {
      // f(x) = op(a x0 + b x0 x1) * x1 + op(x0)^2 : nonlinear arguments exercise the chain rule to second order
      NlProblem p;
      p.number_variables = 2;
      const double middle = 0.5 * (c.low + c.high), half = 0.25 * (c.high - c.low);
      const Expression argument = num(middle) + num(half) * var(0) * var(1);
      p.objectives.push_back({op(c.opcode, {argument}) * var(1) + pow(op(c.opcode, {num(middle) + num(half) * var(0)}), num(2.)), {}});
      p.constraints.push_back({op(c.opcode, {num(middle) + num(half) * var(1)}) / (num(2.) + var(0) * var(0)), {}});
      p.constraint_lower = {-INFINITY};
      p.constraint_upper = {0.};
      const cppasl::NlModel model = testing::write_and_read(p, Format::Text, std::string("unary_") + c.name);
      for (int trial = 0; trial < 3; ++trial) {
         const std::vector<double> x = testing::random_vector(random, 2, -0.9, 0.9);
         const double error = CHECK_DERIVATIVES(model, x, 0.7, std::vector<double>{1.3}, 2e-6);
         if (error > 2e-6) std::fprintf(stderr, "    operation %s\n", c.name);
      }
   }
}

TEST_CASE(binary_and_nary_operations) {
   NlProblem p;
   p.number_variables = 4;
   const Expression x0 = var(0), x1 = var(1), x2 = var(2), x3 = var(3);
   p.objectives.push_back({op_list(54, {
      x0 / (x1 + num(3.)),                      // divide
      pow(x0 + num(2.), x1),                    // power with variable exponent
      pow(num(1.5), x2 * x3),                   // constant base
      pow(x2 + num(2.), num(3.5)),              // constant exponent
      op(48, {x0 + num(0.1), x1 * x2 - num(2.)}), // atan2
      op(80, {x3 - num(0.3), num(2.5)}),        // signpow
      op_list(11, {x0 * x1, x2 * x2, x3}),      // min
      op_list(12, {x0, x1 * x3, num(-5.)}),     // max
      op(35, {op(22, {x0, x1}), x0 * x0 * x2, op(44, {x3})}), // if-then-else
      piecewise_linear({-1., 0.5, 2.}, {-0.5, 0.5}, x0 * x2), // piecewise linear
      op(4, {x0 * num(10.), num(3.)}) * x1,      // remainder (zero second derivatives almost everywhere)
      op(13, {x3 * num(5.)}) * x2 * x2,          // floor
   }), {{1, 0.25}}});
   p.constraints.push_back({op(0, {x0 * x1 * x2 * x3, op(43, {x0 * x0 + num(1.)})}), {}});
   p.constraints.push_back({op(1, {op(46, {x2 * x3}), op(44, {x0 * x1})}), {{3, 2.}}});
   p.constraint_lower = {-INFINITY, 0.};
   p.constraint_upper = {10., 0.};
   const cppasl::NlModel model = testing::write_and_read(p, Format::Binary, "binary_nary");
   std::mt19937_64 random(12);
   for (int trial = 0; trial < 20; ++trial) {
      const std::vector<double> x = testing::random_vector(random, 4, -0.9, 0.9);
      CHECK_DERIVATIVES(model, x, 1.1, std::vector<double>({-0.4, 2.2}), 1e-5);
   }
   // spot-check a value
   cppasl::EvaluationWorkspace workspace(model);
   const std::vector<double> x{0.2, 0.4, 0.3, -0.1};
   double expected = 0.2 / 3.4 + std::pow(2.2, 0.4) + std::pow(1.5, 0.3 * -0.1) + std::pow(2.3, 3.5) + std::atan2(0.3, 0.4 * 0.3 - 2.) +
      std::copysign(std::pow(0.4, 2.5), -0.4) + std::min({0.08, 0.09, -0.1}) + std::max({0.2, -0.04, -5.}) + 0.2 * 0.2 * 0.3;
   const double z = 0.2 * 0.3; // piecewise linear: slope 0.5 on (-0.5, 0.5]
   expected += 0.5 * z + std::fmod(2., 3.) * 0.4 + std::floor(-0.5) * 0.09 + 0.25 * 0.4;
   CHECK_CLOSE(model.evaluate_objective(workspace, x.data()), expected, 1e-14);
}

TEST_CASE(defined_variables_shared_by_functions) {
   NlProblem p;
   p.number_variables = 3;
   p.defined_variables.push_back({{{0, 1.}, {2, -2.}}, op(41, {var(1)})});   // v3 = x0 - 2 x2 + sin(x1)
   p.defined_variables.push_back({{}, var(3) * var(3) + op(44, {var(3)})}); // v4 = v3^2 + exp(v3)
   p.objectives.push_back({var(4) * var(0), {}});
   p.constraints.push_back({var(4) + var(3), {}});
   p.constraints.push_back({op(43, {num(2.) + var(4)}), {{1, 1.}}});
   p.constraint_lower = {-INFINITY, -INFINITY};
   p.constraint_upper = {1., 1.};
   const cppasl::NlModel model = testing::write_and_read(p, Format::Text, "defined_variables");
   std::mt19937_64 random(13);
   const std::vector<double> x = testing::random_vector(random, 3);
   cppasl::EvaluationWorkspace workspace(model);
   const double v3 = x[0] - 2. * x[2] + std::sin(x[1]), v4 = v3 * v3 + std::exp(v3);
   CHECK_CLOSE(model.evaluate_objective(workspace, x.data()), v4 * x[0], 1e-14);
   std::vector<double> c(2);
   model.evaluate_constraints(workspace, x.data(), c.data());
   CHECK_CLOSE(c[0], v4 + v3, 1e-14);
   CHECK_CLOSE(c[1], std::log(2. + v4) + x[1], 1e-14);
   CHECK_DERIVATIVES(model, x, 0.9, std::vector<double>({1.7, -0.6}), 1e-6);
}

TEST_CASE(workspace_caching_and_threads) {
   // alternating points and workspaces must never return stale values
   NlProblem p;
   p.number_variables = 2;
   p.objectives.push_back({op(44, {var(0) * var(1)}), {}});
   p.constraints.push_back({op(41, {var(0)}) * var(1), {}});
   p.constraint_lower = {-INFINITY};
   p.constraint_upper = {0.};
   const cppasl::NlModel model = testing::write_and_read(p, Format::Text, "caching");
   cppasl::EvaluationWorkspace a(model), b(model);
   std::vector<double> x{0.3, 0.7}, y{-0.2, 0.5}, g(2), c(1);
   for (int k = 0; k < 4; ++k) {
      const std::vector<double>& point = (k % 2 == 0) ? x : y;
      cppasl::EvaluationWorkspace& workspace = (k < 2) ? a : b;
      CHECK_CLOSE(model.evaluate_objective(workspace, point.data()), std::exp(point[0] * point[1]), 1e-15);
      model.evaluate_objective_gradient(workspace, point.data(), g.data());
      CHECK_CLOSE(g[0], point[1] * std::exp(point[0] * point[1]), 1e-15);
      model.evaluate_constraints(workspace, point.data(), c.data());
      CHECK_CLOSE(c[0], std::sin(point[0]) * point[1], 1e-15);
   }
   // in-place modification of the same buffer
   cppasl::EvaluationWorkspace w(model);
   std::vector<double> z{0.1, 0.2};
   (void)model.evaluate_objective(w, z.data());
   z[1] = 0.9;
   CHECK_CLOSE(model.evaluate_objective(w, z.data()), std::exp(0.09), 1e-15);
}

TEST_CASE(group_elements) {
   // phi(sum_k t_k) with more than 8 variables uses the group formula phi'' g g^T + phi' sum H_k;
   // operands share variables, contain products and constants, and the chains mix unary and constant-binary steps
   NlProblem p;
   p.number_variables = 14;
   std::vector<Expression> sines{num(1.)}, exponentials, mixed;
   for (int i = 0; i < 14; ++i) sines.push_back(pow(op(41, {var(i)}), num(2.)));
   for (int i = 0; i < 12; ++i) exponentials.push_back(op(44, {var(i) * var((i + 3) % 14) * num(0.5)}));
   for (int i = 0; i < 10; ++i) mixed.push_back(var(i) * var(i + 1) + op(46, {var(i + 4)}));
   mixed.push_back(num(2.));
   p.objectives.push_back({op_list(54, {op(39, {op_list(54, sines)}),                          // sqrt(1 + sum sin^2)
      op(43, {op_list(54, exponentials)}),                                                         // log-sum-exp
      num(0.3) * pow(num(2.) - op_list(54, mixed) / num(7.), num(3.)) }), {}});                   // (2 - sum/7)^3
   p.constraints.push_back({op(37, {num(0.1) * op_list(54, exponentials)}), {}});
   p.constraints.push_back({op(44, {neg(op_list(54, sines))}), {{3, 1.}}});
   p.constraint_lower = {-INFINITY, -INFINITY};
   p.constraint_upper = {1., 1.};
   const cppasl::NlModel model = testing::write_and_read(p, Format::Text, "groups");
   cppasl::ReaderOptions no_detection;
   std::mt19937_64 random(15);
   for (int trial = 0; trial < 5; ++trial) {
      const std::vector<double> x = testing::random_vector(random, 14, -0.8, 0.8);
      CHECK_DERIVATIVES(model, x, 1.2, std::vector<double>({0.8, -1.5}), 1e-6);
   }
}

TEST_CASE(shared_defined_variables_chain_rule) {
   // defined variables evaluated once and shared by several functions: nested ones, one inside a piecewise-constant
   // operation (no adjoint), one in a group element, a second objective outside the Lagrangian
   NlProblem p;
   p.number_variables = 12;
   std::vector<Expression> exponentials;
   for (int j = 0; j < 12; ++j) exponentials.push_back(op(44, {var(j) * num(0.3)}));
   p.defined_variables.push_back({{{0, 1.}}, op(43, {op_list(54, exponentials)})});        // v12 = x0 + log(sum exp)
   p.defined_variables.push_back({{{3, -2.}}, op(41, {var(12) * var(1)}) + var(2) * var(2)}); // v13 = -2 x3 + sin(v12 x1) + x2^2
   std::vector<Expression> group_terms;
   for (int j = 0; j < 10; ++j) group_terms.push_back(pow(var(j) - var(13) * num(0.1), num(2.)));
   p.objectives.push_back({op(39, {num(1.) + op_list(54, group_terms)}) + var(12) * var(13), {{5, 1.}}});
   p.objectives.push_back({op(44, {var(13)}), {}});
   p.objective_senses = {0, 0};
   p.constraints.push_back({var(12) * var(4) + op(46, {var(13)}), {}});
   p.constraints.push_back({op(13, {var(12)}) * var(6) + pow(var(13), num(3.)), {{7, 1.}}});
   p.constraints.push_back({op(43, {num(5.) + var(12)}) * var(13) * var(11), {}});
   p.constraint_lower = {-INFINITY, -INFINITY, -INFINITY};
   p.constraint_upper = {1., 1., 1.};
   const cppasl::NlModel model = testing::write_and_read(p, Format::Text, "shared_defined");
   std::mt19937_64 random(17);
   for (int trial = 0; trial < 5; ++trial) {
      const std::vector<double> x = testing::random_vector(random, 12, -0.7, 0.7);
      CHECK_DERIVATIVES(model, x, 0.9, std::vector<double>({1.1, -0.7, 0.4}), 2e-6);
   }
}
