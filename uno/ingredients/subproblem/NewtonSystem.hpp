// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NEWTONSYSTEM_H
#define UNO_NEWTONSYSTEM_H

#include <cstddef>
#include "NewtonSystemContext.hpp"
#include "Subproblem.hpp"
#include "linear_algebra/View.hpp"
#include "../interfaces/C/uno_int.h"

namespace uno {
   // Newton system of the equality-constrained subproblem. The optimization problems declare the uncondensed (possibly
   // unsymmetric) system as a symbolic block system (see OptimizationProblem::create_newton_system and
   // SymbolicNewtonSystem.hpp), with blocks of unknowns (Δx, -Δy, ...) and a list of blocks to condense. The subproblem
   // solvers see the condensed symmetric system
   //    [H + δ_w I    *    ]
   //    [J          -δ_c I ]
   // where H is the curvature block (+ the Schur complements of the condensed blocks), J the constraint Jacobian and δ_w,
   // δ_c the regularization blocks set by the inertia correction strategy.
   // Call protocol, at each new iterate: evaluate() -> assemble_rhs() -> back_substitute() (possibly several times, e.g.
   // for second-order corrections). It is checked: a violation throws std::logic_error
   class NewtonSystem {
   public:
      virtual ~NewtonSystem() = default;

      // condensed system
      [[nodiscard]] virtual size_t dimension() const = 0;
      [[nodiscard]] virtual size_t number_nonzeros() const = 0;
      // sparsity pattern (lower triangle) of the symmetric condensed matrix
      virtual void compute_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, uno_int solver_indexing) const = 0;
      // values of the condensed matrix at the current iterate (the regularization blocks are left untouched)
      virtual void evaluate(const NewtonSystemContext& context, View<double> values) const = 0;
      // views on the blocks of the matrix values (Hessian block, primal regularization, Jacobian, dual regularization)
      [[nodiscard]] virtual BlockAugmentedMatrix block_augmented_matrix(View<double> values) const = 0;
      // right-hand side of the condensed system at the current iterate (precondition: evaluate() was called at this iterate)
      virtual void assemble_rhs(const NewtonSystemContext& context, View<double> rhs) const = 0;
      // eliminated unknowns (number_eliminated() of them, ordered by condensed block) from the solution of the condensed
      // system (precondition: assemble_rhs() was called at the current iterate)
      [[nodiscard]] virtual size_t number_eliminated() const = 0;
      virtual void back_substitute(View<const double> solution, View<double> eliminated_solution) const = 0;
   };
} // namespace

#endif // UNO_NEWTONSYSTEM_H
