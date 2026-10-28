// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NEWTONSYSTEM_H
#define UNO_NEWTONSYSTEM_H

#include <cstddef>
#include <utility>
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

      // uncondensed system (iterative refinement), with vectors laid out as [condensed unknowns | eliminated unknowns].
      // Preconditions: evaluate() (and assemble_rhs() for uncondense_rhs) at the current iterate
      [[nodiscard]] virtual size_t uncondensed_dimension() const = 0;
      // nonzeros [begin, end) of the Schur complement in the values of the condensed matrix: the product of the condensed
      // matrix without them is the product of the kept block of the uncondensed matrix (regularization included)
      [[nodiscard]] virtual std::pair<size_t, size_t> schur_complement_nonzeros() const = 0;
      // uncondensed right-hand side of a condensed right-hand side whose eliminated rows are those of assemble_rhs()
      virtual void uncondense_rhs(View<const double> condensed_rhs, View<double> uncondensed_rhs) const = 0;
      virtual void condense_rhs(View<const double> uncondensed_rhs, View<double> condensed_rhs) const = 0;
      // uncondensed solution, given an uncondensed right-hand side and the solution of the condensed system
      virtual void expand_solution(View<const double> uncondensed_rhs, View<const double> condensed_solution,
         View<double> uncondensed_solution) const = 0;
      // on entry, the condensed entries of product contain the product of the kept block (condensed matrix without the
      // Schur complement) with the condensed entries of vector; on exit, product is the uncondensed product
      virtual void complete_uncondensed_product(View<const double> vector, View<double> product) const = 0;
   };
} // namespace

#endif // UNO_NEWTONSYSTEM_H
