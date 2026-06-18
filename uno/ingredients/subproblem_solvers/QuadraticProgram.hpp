// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_QUADRATICPROGRAM_H
#define UNO_QUADRATICPROGRAM_H

#include <cstddef>
#include <vector>
#include "../interfaces/C/uno_int.h"

namespace uno {
   // forward declarations
   class Evaluations;
   class SolverWorkspace;
   class Statistics;
   class Subproblem;
   template <typename ElementType>
   class Vector;
   class WarmstartInformation;

   // Subproblem-free description of the LP/QP solved by an LPSolver/QPSolver backend.
   //
   // It is abstract because each backend stores the data in its own native sparse format (BQPD's
   // packed gradients + weak-CSR Jacobian, HiGHS's CSC HighsModel). There are two ways to populate it:
   //   - from a Subproblem, via build(Statistics&, const Subproblem&, ...): IQPSolver uses this in the
   //     full solver, after allocating storage with initialize_memory(const Subproblem&);
   //   - directly from data, via the COO build() overload below: the objective gradient is dense and
   //     the constraint Jacobian and Lagrangian Hessian are passed in COO (coordinate) format, which is
   //     the uniform exchange format. The QuadraticProgram converts COO to the backend's preferred
   //     layout internally (weak-CSR for BQPD, CSC for HiGHS). This path needs no Subproblem and is what
   //     the functional tests use; it infers the dimensions from the data and allocates as needed.
   //
   // Either way, the backend's solve() reads only the QuadraticProgram, never a Subproblem.
   class QuadraticProgram {
   public:
      QuadraticProgram() = default;
      virtual ~QuadraticProgram();

      // dimensions, set by whichever population path runs (initialize_memory or the COO build)
      size_t number_variables{0};
      size_t number_constraints{0};
      size_t number_jacobian_nonzeros{0};

      // allocate native storage and compute the (iteration-invariant) sparsity patterns from a Subproblem
      virtual void initialize_memory(const Subproblem& subproblem) = 0;

      // populate from a Subproblem at the current iterate (storage allocated by initialize_memory)
      virtual void fill(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) = 0;

      // populate directly from data: dense objective gradient, COO constraint Jacobian (row = constraint,
      // column = variable), COO Lagrangian Hessian (one triangle; empty for an LP), and concatenable
      // variable/constraint bounds. Dimensions are inferred and native storage is (re)allocated.
      virtual void fill(const Vector<double>& linear_objective,
         const Vector<uno_int>& jacobian_row_indices, const Vector<uno_int>& jacobian_column_indices,
         const Vector<double>& jacobian_values,
         const Vector<uno_int>& hessian_row_indices, const Vector<uno_int>& hessian_column_indices,
         const Vector<double>& hessian_values,
         const std::vector<double>& variables_lower_bounds, const std::vector<double>& variables_upper_bounds,
         const std::vector<double>& constraints_lower_bounds, const std::vector<double>& constraints_upper_bounds) = 0;

      // bridge to the workspace used by the caller for predicted-reduction quadratic forms
      [[nodiscard]] virtual SolverWorkspace& get_workspace() = 0;
   };
} // namespace

#endif // UNO_QUADRATICPROGRAM_H
