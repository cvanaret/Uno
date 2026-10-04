#pragma once

#include <string>

namespace cppasl {

   /// Structural class of a problem, detected exactly from the expression graphs.
   enum class ProblemType {
      LinearProgram,                            ///< linear objective and constraints
      QuadraticProgram,                         ///< quadratic objective, linear constraints
      QuadraticallyConstrainedQuadraticProgram, ///< >= 1 quadratic constraint, no other nonlinearity
      NonlinearProgram                          ///< >= 1 non-quadratic nonlinear function
   };

   enum class Convexity { Convex, Nonconvex, Unknown };

   struct ClassificationOptions {
      /// Quadratic forms supported on more variables than this are not factorized: their convexity is
      /// Unknown unless a cheap certificate (negative 1x1/2x2 minor, diagonal dominance) applies.
      int maximum_dense_dimension{4000};
      /// H is declared positive semidefinite iff H + tolerance * max(1, max_i |H_ii|) I is Cholesky-factorizable.
      double semidefiniteness_tolerance{1e-8};
   };

   struct ProblemClassification {
      ProblemType type{ProblemType::LinearProgram};
      Convexity convexity{Convexity::Unknown};
      bool has_integer_variables{false};

      /// e.g. "LP", "convex QP", "nonconvex QCQP", "MINLP"
      [[nodiscard]] std::string to_string() const;
   };

} // namespace cppasl
