// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_QUADRATICPROGRAM_H
#define UNO_QUADRATICPROGRAM_H

#include <cstddef>

namespace uno {
   // forward declarations
   class Evaluations;
   class SolverWorkspace;
   class Statistics;
   class Subproblem;
   class WarmstartInformation;

   // Subproblem-free description of the LP/QP solved by an LPSolver/QPSolver backend.
   //
   // It is abstract because each backend stores the data in its own native sparse format (BQPD's
   // packed gradients + weak-CSR Jacobian, HiGHS's CSC HighsModel). IQPSolver owns the backend,
   // retrieves the backend's QuadraticProgram via LPSolver::get_quadratic_program(), and populates it
   // from a Subproblem with build(). The backend's solve() then reads only the QuadraticProgram, never
   // a Subproblem, which is what makes the backends testable on hand-built data.
   //
   // An LP is the curvature-free case: build() simply leaves the Hessian empty. The capability
   // asymmetry (a QP solver can solve LPs but not vice versa) lives in the QPSolver : LPSolver
   // hierarchy, not here.
   class QuadraticProgram {
   public:
      QuadraticProgram(size_t number_variables, size_t number_constraints);
      virtual ~QuadraticProgram();

      const size_t number_variables;
      const size_t number_constraints;

      // allocate native storage and compute the (iteration-invariant) sparsity patterns
      virtual void initialize_memory(const Subproblem& subproblem) = 0;

      // evaluate the objective gradient, constraint Jacobian, bounds and Hessian at the current iterate.
      // warmstart_information indicates which parts changed since the previous build.
      virtual void build(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) = 0;

      // bridge to the workspace used by the caller for predicted-reduction quadratic forms
      [[nodiscard]] virtual SolverWorkspace& get_workspace() = 0;
   };
} // namespace

#endif // UNO_QUADRATICPROGRAM_H
