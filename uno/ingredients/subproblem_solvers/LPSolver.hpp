// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LPSOLVER_H
#define UNO_LPSOLVER_H

namespace uno {
   // forward declarations
   class Direction;
   class Statistics;
   class Subproblem;
   template <typename ElementType>
   class Vector;
   class WarmstartInformation;

   class LPSolver {
   public:
      LPSolver() = default;
      virtual ~LPSolver() = default;

      virtual void initialize_memory(const Subproblem& subproblem) = 0;

      virtual void solve(Statistics& statistics, Subproblem& subproblem, const Vector<double>& initial_point,
         Direction& direction, const WarmstartInformation& warmstart_information) = 0;

      virtual void compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const = 0;
      virtual void compute_constraint_jacobian_transposed_vector_product(const Vector<double>& vector, Vector<double>& result) const = 0;
      [[nodiscard]] virtual double compute_hessian_quadratic_product(const Vector<double>& vector) const = 0;
   };
} // namespace

#endif // UNO_LPSOLVER_H