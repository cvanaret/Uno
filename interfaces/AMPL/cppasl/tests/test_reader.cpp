// Reader tests: hand-written text file, text vs binary encodings, bounds, suffixes, integrality, defined variables.
#include <cstring>
#include <fstream>
#include <limits>
#include "test_utils.hpp"

using namespace nlwriter;

namespace {
   // HS071 as AMPL writes it in text mode (with comments), plus a linear term x2 inside the objective expression
   const char* hs071 =
      "g3 1 1 0\t# problem hs071\n"
      " 4 2 1 0 1\t# vars, constraints, objectives, ranges, eqns\n"
      " 2 1\t# nonlinear constraints, objectives\n"
      " 0 0\t# network constraints: nonlinear, linear\n"
      " 4 4 4\t# nonlinear vars in constraints, objectives, both\n"
      " 0 0 0 1\t# linear network variables; functions; arith, flags\n"
      " 0 0 0 0 0\t# discrete variables: binary, integer, nonlinear (b,c,o)\n"
      " 8 4\t# nonzeros in Jacobian, gradients\n"
      " 0 0\t# max name lengths: constraints, variables\n"
      " 0 0 0 0 0\t# common exprs: b,c,o,c1,o1\n"
      "C0\t#c1\no2\t#*\no2\no2\nv0\nv1\nv2\nv3\n"
      "C1\t#c2\no54\t#sumlist\n4\no5\t#^\nv0\nn2\no5\nv1\nn2\no5\nv2\nn2\no5\nv3\nn2\n"
      "O0 0\t#obj\no0\t#+\no2\no2\nv0\nv3\no54\n3\nv0\nv1\nv2\nv2\n"
      "x4\t# initial guess\n0 1\n1 5\n2 5\n3 1\n"
      "r\t#2 ranges (rhs's)\n2 25\n4 40\n"
      "b\t#4 bounds (on variables)\n0 1 5\n0 1 5\n0 1 5\n0 1 5\n"
      "k3\t#intermediate Jacobian column lengths\n2\n4\n6\n"
      "J0 4\n0 0\n1 0\n2 0\n3 0\n"
      "J1 4\n0 0\n1 0\n2 0\n3 0\n"
      "G0 4\n0 0\n1 0\n2 0\n3 0\n";
} // namespace

TEST_CASE(reads_handwritten_hs071) {
   const cppasl::NlModel model = cppasl::NlModel::read_from_memory(hs071, std::strlen(hs071));
   CHECK(model.number_variables() == 4);
   CHECK(model.number_constraints() == 2);
   CHECK(model.header().format == cppasl::NlFormat::Text);
   CHECK(model.initial_primal_point() == std::vector<double>({1., 5., 5., 1.}));
   CHECK(model.constraint_lower_bounds()[0] == 25. && std::isinf(model.constraint_upper_bounds()[0]));
   CHECK(model.constraint_lower_bounds()[1] == 40. && model.constraint_upper_bounds()[1] == 40.);
   CHECK(model.variable_lower_bounds()[2] == 1. && model.variable_upper_bounds()[2] == 5.);
   CHECK(model.number_jacobian_nonzeros() == 8);
   cppasl::EvaluationWorkspace workspace(model);
   const std::vector<double>& x = model.initial_primal_point();
   CHECK_CLOSE(model.evaluate_objective(workspace, x.data()), 16., 1e-15);
   std::vector<double> gradient(4), constraints(2);
   model.evaluate_objective_gradient(workspace, x.data(), gradient.data());
   CHECK(gradient == std::vector<double>({12., 1., 2., 11.}));
   model.evaluate_constraints(workspace, x.data(), constraints.data());
   CHECK(constraints == std::vector<double>({25., 52.}));
   CHECK(model.constraint_structure(0) == cppasl::FunctionStructure::Nonlinear);
   CHECK(model.constraint_structure(1) == cppasl::FunctionStructure::Quadratic);
   CHECK(model.objective_structure() == cppasl::FunctionStructure::Nonlinear);
   CHECK(model.classify().to_string() == "NLP");
   CHECK(model.number_hessian_nonzeros() == 10); // dense lower triangle of a 4x4 matrix
   CHECK_DERIVATIVES(model, x, 1.3, std::vector<double>({0.7, -2.1}), 1e-6);
}

TEST_CASE(rejects_malformed_files) {
   auto fails = [](const std::string& contents) {
      try {
         cppasl::NlModel::read_from_memory(contents.data(), contents.size());
      }
      catch (const std::runtime_error&) {
         return true;
      }
      return false;
   };
   std::string truncated(hs071);
   truncated.resize(truncated.find("C1"));
   truncated += "C1\no54\n4\nv0\n"; // incomplete sum
   CHECK(fails(truncated));
   std::string bad_variable(hs071);
   bad_variable.replace(bad_variable.find("v3\nC1"), 2, "v9");
   CHECK(fails(bad_variable));
   CHECK(fails("x3 1 1 0\n"));
}

namespace {
   /// random NLP touching most features of the format
   NlProblem feature_problem() {
      NlProblem p;
      p.number_variables = 6;
      p.defined_variables.push_back({{{0, 2.}, {1, -1.}}, op(44, {var(2) * num(0.5)})}); // v6 = 2x0 - x1 + exp(x2/2)
      p.defined_variables.push_back({{}, var(6) * var(3)});                                 // v7 = v6 * x3
      p.objectives.push_back({op(43, {num(3.) + pow(var(4), num(2.))}) + var(7), {{5, 1.5}, {0, -1.}}});
      p.objective_senses = {0};
      p.constraints.push_back({op(41, {var(0) * var(1)}) + op(46, {var(7)}), {{2, 1.}}});
      p.constraints.push_back({op(39, {num(1.) + var(3) * var(3)}), {{4, 3.}}});
      p.constraints.push_back({{}, {{0, 1.}, {1, 1.}, {5, -2.}}});
      p.constraints.push_back({{}, {{3, 4.}}});
      p.constraint_lower = {-1., -INFINITY, 0., 1.};
      p.constraint_upper = {1., 10., INFINITY, 1.};
      p.variable_lower = {-1., -INFINITY, 0., -5., -INFINITY, 2.};
      p.variable_upper = {1., 3., INFINITY, 5., INFINITY, 2.};
      p.initial_point = {{0, 0.5}, {3, -0.25}, {5, 2.}};
      p.initial_dual = {{1, 0.125}};
      p.suffixes.push_back({0, "priority", {{1, 3.}, {4, 7.}}});
      p.suffixes.push_back({5, "scale", {{2, 0.5}}});
      return p;
   }
} // namespace

TEST_CASE(text_and_binary_encodings_agree) {
   const NlProblem problem = feature_problem();
   const cppasl::NlModel text = testing::write_and_read(problem, Format::Text, "features");
   const cppasl::NlModel binary = testing::write_and_read(problem, Format::Binary, "features");
   CHECK(binary.header().format == cppasl::NlFormat::Binary);
   CHECK(text.variable_lower_bounds() == binary.variable_lower_bounds());
   CHECK(text.variable_upper_bounds() == binary.variable_upper_bounds());
   CHECK(text.constraint_lower_bounds() == binary.constraint_lower_bounds());
   CHECK(text.constraint_upper_bounds() == binary.constraint_upper_bounds());
   CHECK(text.initial_primal_point() == binary.initial_primal_point());
   CHECK(text.initial_dual_point() == binary.initial_dual_point());
   CHECK(text.jacobian_column_indices() == binary.jacobian_column_indices());
   CHECK(text.hessian_row_indices() == binary.hessian_row_indices());
   CHECK(text.hessian_column_starts() == binary.hessian_column_starts());
   CHECK(text.initial_primal_point() == std::vector<double>({0.5, 0., 0., -0.25, 0., 2.}));
   CHECK(text.initial_dual_point()[1] == 0.125);
   CHECK(std::isinf(text.variable_lower_bounds()[1]) && text.variable_lower_bounds()[5] == 2.);
   // suffixes
   for (const cppasl::NlModel* model: {&text, &binary}) {
      CHECK(model->suffixes().size() == 2);
      CHECK(model->suffixes()[0].name == "priority" && model->suffixes()[0].values[4] == 7. && !model->suffixes()[0].is_real);
      CHECK(model->suffixes()[1].target == cppasl::SuffixTarget::Constraints && model->suffixes()[1].is_real);
      CHECK(model->suffixes()[1].values[2] == 0.5);
   }
   // identical evaluations (bitwise: same shortest round-trip decimal representation)
   std::mt19937_64 random(1);
   const std::vector<double> x = testing::random_vector(random, 6, 0.1, 0.9), y = testing::random_vector(random, 4);
   cppasl::EvaluationWorkspace text_workspace(text), binary_workspace(binary);
   CHECK(text.evaluate_objective(text_workspace, x.data()) == binary.evaluate_objective(binary_workspace, x.data()));
   std::vector<double> ct(4), cb(4), ht(text.number_hessian_nonzeros()), hb(binary.number_hessian_nonzeros());
   text.evaluate_constraints(text_workspace, x.data(), ct.data());
   binary.evaluate_constraints(binary_workspace, x.data(), cb.data());
   CHECK(ct == cb);
   text.evaluate_lagrangian_hessian(text_workspace, x.data(), 2., y.data(), ht.data());
   binary.evaluate_lagrangian_hessian(binary_workspace, x.data(), 2., y.data(), hb.data());
   CHECK(ht == hb);
   CHECK_DERIVATIVES(text, x, 2., y, 1e-6);
   CHECK(text.constraint_structure(2) == cppasl::FunctionStructure::Linear);
   CHECK(text.constraint_structure(1) == cppasl::FunctionStructure::Nonlinear);
}

TEST_CASE(integer_variables_and_complementarity) {
   NlProblem p;
   p.number_variables = 5;
   p.number_binary_variables = 1;
   p.number_integer_variables = 2;
   p.objectives.push_back({{}, {{0, 1.}, {4, 2.}}});
   p.objective_senses = {1};
   p.constraints.push_back({{}, {{0, 1.}, {1, 1.}, {2, 1.}}});
   p.constraint_lower = {0.};
   p.constraint_upper = {4.};
   std::string contents = p.write(Format::Text);
   // turn the constraint bounds into a complementarity with variable 2 (1-based: 2 -> index 1)
   contents.replace(contents.find("\nr\n0 0 4\n"), 9, "\nr\n5 1 2\n");
   const cppasl::NlModel model = cppasl::NlModel::read_from_memory(contents.data(), contents.size());
   CHECK(!model.is_integer_variable(0) && !model.is_integer_variable(1));
   CHECK(model.is_integer_variable(2) && model.is_integer_variable(3) && model.is_integer_variable(4));
   CHECK(model.is_maximization());
   CHECK(model.complementary_variables()[0] == 1);
   CHECK(model.constraint_lower_bounds()[0] == 0. && std::isinf(model.constraint_upper_bounds()[0]));
   CHECK(model.classify().to_string() == "convex MILP");
}

TEST_CASE(writes_ampl_solution_file) {
   const cppasl::NlModel model = cppasl::NlModel::read_from_memory(hs071, std::strlen(hs071));
   const std::vector<double> x{1., 4.743, 3.821, 1.379}, y{-0.552, 0.161};
   const std::string file_name = "/tmp/cppasl_test_hs071.sol";
   model.write_solution(file_name, "cppasl: optimal solution\nsecond line", x.data(), y.data(), 0);
   std::ifstream file(file_name);
   std::string contents((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
   CHECK(contents == "cppasl: optimal solution\nsecond line\n\nOptions\n3\n1\n1\n0\n2\n2\n4\n4\n"
      "-0.552\n0.161\n1\n4.743\n3.821\n1.379\nobjno 0 0\n");
}
