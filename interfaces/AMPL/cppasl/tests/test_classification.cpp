// LP / QP / QCQP / NLP classification and convexity detection.
#include "test_utils.hpp"

using namespace nlwriter;

namespace {
   std::string classify(const Expression& objective, int sense, const std::vector<Expression>& constraints,
         const std::vector<std::pair<double, double>>& bounds, cppasl::ClassificationOptions options = {}) {
      NlProblem p;
      p.number_variables = 4;
      p.objectives.push_back({objective, {{3, 1.}}});
      p.objective_senses = {sense};
      for (std::size_t i = 0; i < constraints.size(); ++i) {
         p.constraints.push_back({constraints[i], {{0, 1.}}});
         p.constraint_lower.push_back(bounds[i].first);
         p.constraint_upper.push_back(bounds[i].second);
      }
      return testing::write_and_read(p, Format::Text, "classification").classify(options).to_string();
   }

   const Expression x0 = var(0), x1 = var(1), x2 = var(2);
   const Expression convex = x0 * x0 + x0 * x1 + x1 * x1 + num(2.) * x2 * x2;   // [[2,1],[1,2]] and 4: PD
   const Expression semidefinite = pow(x0 - x1, num(2.)) + x2 * x2;               // PSD, singular
   const Expression indefinite = x0 * x0 - x1 * x1;
   const Expression hidden_indefinite = op_list(54, {x0 * x0, x1 * x1, x2 * x2, num(1.9) * x0 * x1, num(1.9) * x1 * x2,
      num(-1.9) * x0 * x2}); // every 2x2 minor positive, but not PSD (x = (1,-1,1)): needs the factorization
} // namespace

TEST_CASE(classification_of_problem_types) {
   const double inf = INFINITY;
   CHECK(classify({}, 0, {{}}, {{0., 1.}}) == "convex LP");
   CHECK(classify(convex, 0, {{}}, {{0., 1.}}) == "convex QP");
   CHECK(classify(semidefinite, 0, {}, {}) == "convex QP");
   CHECK(classify(indefinite, 0, {}, {}) == "nonconvex QP");
   CHECK(classify(hidden_indefinite, 0, {}, {}) == "nonconvex QP");
   CHECK(classify(neg(convex), 1, {}, {}) == "convex QP");      // maximize a concave quadratic
   CHECK(classify(convex, 1, {}, {}) == "nonconvex QP");        // maximize a convex quadratic
   CHECK(classify(convex, 0, {convex}, {{-inf, 1.}}) == "convex QCQP");
   CHECK(classify(convex, 0, {neg(convex)}, {{-1., inf}}) == "convex QCQP");
   CHECK(classify(convex, 0, {convex}, {{1., inf}}) == "nonconvex QCQP");
   CHECK(classify(convex, 0, {convex}, {{1., 1.}}) == "nonconvex QCQP");   // quadratic equality
   CHECK(classify(convex, 0, {x0 * x1}, {{-inf, 1.}}) == "nonconvex QCQP"); // bilinear
   CHECK(classify({}, 0, {semidefinite, convex}, {{-inf, 1.}, {-inf, 2.}}) == "convex QCQP");
   CHECK(classify(convex, 0, {op(44, {x0})}, {{-inf, 1.}}) == "NLP");
   CHECK(classify(op(43, {x0}), 0, {}, {}) == "NLP");
   // the dense factorization is skipped above the size limit
   cppasl::ClassificationOptions small;
   small.maximum_dense_dimension = 2;
   CHECK(classify(hidden_indefinite, 0, {}, {}, small) == "QP");
}

TEST_CASE(classification_of_large_block_structured_forms) {
   // sum of 20000 independent 2x2 PSD blocks (x_2k - x_2k+1)^2: components keep the factorizations tiny;
   // flipping one block to indefinite is detected
   for (bool flip: {false, true}) {
      NlProblem p;
      p.number_variables = 40000;
      std::vector<Expression> terms;
      for (int k = 0; k < 20000; ++k) {
         const Expression block = num(4.) * var(2 * k) * var(2 * k) + num(4.) * var(2 * k) * var(2 * k + 1) +
            num(1.) * var(2 * k + 1) * var(2 * k + 1) + num(1.) * var(2 * k + 1) * var(2 * k + 1); // [[8,4],[4,4]]
         terms.push_back((flip && k == 12345) ? neg(block) : block);
      }
      p.objectives.push_back({op_list(54, terms), {}});
      const cppasl::NlModel model = testing::write_and_read(p, Format::Binary, "blocks");
      CHECK(model.classify().to_string() == (flip ? "nonconvex QP" : "convex QP"));
   }
}
