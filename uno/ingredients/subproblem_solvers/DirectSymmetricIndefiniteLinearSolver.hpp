// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DIRECTSYMMETRICINDEFINITELINEARSOLVER_H
#define UNO_DIRECTSYMMETRICINDEFINITELINEARSOLVER_H

#include <cstddef>
#include <memory>
#include <optional>
#include <type_traits>
#include <utility>
#include "SymmetricIndefiniteLinearSolver.hpp"
#include "COOLinearSystem.hpp"
#include "ScalingStrategy.hpp"
#include "ingredients/inertia_correction_strategies/Inertia.hpp"

namespace uno {
   // the public factorization and solves are non-virtual: they apply the optional linear system scaling (S K S y = S b,
   // x = S y) around the solver-specific protected primitives factorize() and solve()
   template <typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver: public SymmetricIndefiniteLinearSolver<ElementType> {
      static_assert(std::is_same_v<ElementType, double>, "The direct linear solvers operate on double");

   public:
      DirectSymmetricIndefiniteLinearSolver() = default;
      ~DirectSymmetricIndefiniteLinearSolver() override = default;

      virtual void do_symbolic_analysis() = 0;

      void do_numerical_factorization(bool is_matrix_positive_definite) {
         if (this->scaling_strategy == nullptr) {
            this->factorize(is_matrix_positive_definite);
            return;
         }
         COOLinearSystem& linear_system = this->get_coo_linear_system();
         const bool scaled = this->scaling_strategy->scale_matrix(linear_system);
         try {
            this->factorize(is_matrix_positive_definite);
         }
         catch (...) {
            if (scaled) {
               this->scaling_strategy->unscale_matrix(linear_system);
            }
            throw;
         }
         // the callers see the unscaled matrix
         if (scaled) {
            this->scaling_strategy->unscale_matrix(linear_system);
         }
         this->scaling_strategy->notify_factorization(this->number_delayed_pivots(), linear_system.dimension);
      }

      void solve_indefinite_system(ElementType* solution) final {
         this->solve_indefinite_system(this->get_linear_system().rhs.data(), solution, 1);
      }

      void solve_indefinite_system(const ElementType* rhs, ElementType* solution, size_t number_of_rhs) final {
         if (this->scaling_strategy != nullptr && this->scaling_strategy->is_last_factorization_scaled()) {
            const size_t dimension = this->get_linear_system().dimension;
            this->solve(this->scaling_strategy->scale_rhs(rhs, dimension, number_of_rhs), solution, number_of_rhs);
            this->scaling_strategy->scale_solution(solution, dimension, number_of_rhs);
         }
         else {
            this->solve(rhs, solution, number_of_rhs);
         }
      }

      [[nodiscard]] virtual Inertia get_inertia() const = 0;
      [[nodiscard]] virtual size_t number_negative_eigenvalues() const = 0;
      // [[nodiscard]] virtual bool matrix_is_positive_definite() const = 0;
      [[nodiscard]] virtual bool matrix_is_singular() const = 0;
      [[nodiscard]] virtual size_t rank() const = 0;

      // increase the pivoting tolerance (to improve the accuracy of the next factorizations) if the solver supports it.
      // Returns false if the tolerance cannot be increased (not supported or maximum value reached).
      // The new tolerance is used from the next numerical factorization on, and is never decreased
      [[nodiscard]] virtual bool increase_pivoting_tolerance() { return false; }

      // number of pivots delayed by numerical pivoting in the last factorization, if reported by the solver
      [[nodiscard]] virtual std::optional<size_t> number_delayed_pivots() const { return std::nullopt; }

      [[nodiscard]] virtual COOLinearSystem& get_coo_linear_system() = 0;

      // linear system scaling (see ScalingStrategy). Without a strategy (e.g. solvers that only compute inertias), the
      // matrices are never scaled
      void set_scaling_strategy(std::unique_ptr<ScalingStrategy> strategy) { this->scaling_strategy = std::move(strategy); }
      [[nodiscard]] const ScalingStrategy* get_scaling_strategy() const { return this->scaling_strategy.get(); }
      // the next factorization is that of a new matrix (as opposed to a refactorization of the current matrix with a
      // different regularization or pivoting tolerance)
      void begin_new_matrix() {
         if (this->scaling_strategy != nullptr) {
            this->scaling_strategy->begin_new_matrix();
         }
      }
      // switch on the scaling from the next factorization on. Returns false if there is no strategy or it is already on
      [[nodiscard]] bool activate_scaling() {
         return this->scaling_strategy != nullptr && this->scaling_strategy->activate();
      }

   protected:
      std::unique_ptr<ScalingStrategy> scaling_strategy{};

      // solver-specific numerical factorization of get_linear_system().matrix_values
      virtual void factorize(bool is_matrix_positive_definite) = 0;
      // solver-specific solve with `number_of_rhs` right-hand sides stored column-major (the two blocks may alias)
      virtual void solve(const ElementType* rhs, ElementType* solution, size_t number_of_rhs) = 0;
   };
} // namespace

#endif // UNO_DIRECTSYMMETRICINDEFINITELINEARSOLVER_H