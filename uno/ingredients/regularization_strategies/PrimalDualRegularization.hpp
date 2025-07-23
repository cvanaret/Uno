// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALDUALREGULARIZATION_H
#define UNO_PRIMALDUALREGULARIZATION_H

#include <cassert>
#include <string>
#include "RegularizationStrategy.hpp"
#include "UnstableRegularization.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "options/Options.hpp"
#include "symbolic/Collection.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   template <typename ElementType>
   class PrimalDualRegularization: public RegularizationStrategy<ElementType> {
   public:
      explicit PrimalDualRegularization(const Options& options);

      void initialize_memory(const OptimizationProblem& problem, const HessianModel& hessian_model) override;
      void initialize_statistics(Statistics& statistics, const Options& options) override;

      void regularize_hessian(Statistics& statistics, const Subproblem& subproblem, SymmetricMatrix<size_t, ElementType>& hessian,
         const Inertia& expected_inertia) override;
      void regularize_hessian(Statistics& statistics, const Subproblem& subproblem, SymmetricMatrix<size_t, ElementType>& hessian,
         const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<size_t, double>& linear_solver) override;
      void regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         SymmetricMatrix<size_t, ElementType>& augmented_matrix, ElementType dual_regularization_parameter,
         const Inertia& expected_inertia) override;
      void regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         SymmetricMatrix<size_t, ElementType>& augmented_matrix, ElementType dual_regularization_parameter,
         const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<size_t, double>& linear_solver) override;

      [[nodiscard]] bool performs_primal_regularization() const override;
      [[nodiscard]] bool performs_dual_regularization() const override;
      [[nodiscard]] double get_primal_regularization_factor() const override;
      [[nodiscard]] std::string get_name() const override;

   protected:
      const std::string& optional_linear_solver_name;
      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<size_t, double>> optional_linear_solver{};
      size_t number_variables{};
      size_t number_constraints{};
      size_t dimension{};
      size_t number_hessian_nonzeros{};
      size_t number_jacobian_nonzeros{};
      ElementType primal_regularization{0.};
      ElementType dual_regularization{0.};
      ElementType previous_primal_regularization{0.};
      const ElementType regularization_failure_threshold;
      const ElementType primal_regularization_initial_factor;
      const ElementType dual_regularization_fraction;
      const ElementType primal_regularization_lb;
      const ElementType primal_regularization_decrease_factor;
      const ElementType primal_regularization_fast_increase_factor;
      const ElementType primal_regularization_slow_increase_factor;
      const size_t threshold_unsuccessful_attempts;
      bool symbolic_analysis_performed{false};
   };

   template <typename ElementType>
   PrimalDualRegularization<ElementType>::PrimalDualRegularization(const Options& options):
         RegularizationStrategy<ElementType>(),
         optional_linear_solver_name(options.get_string("linear_solver")),
         regularization_failure_threshold(ElementType(options.get_double("regularization_failure_threshold"))),
         primal_regularization_initial_factor(ElementType(options.get_double("primal_regularization_initial_factor"))),
         dual_regularization_fraction(ElementType(options.get_double("dual_regularization_fraction"))),
         primal_regularization_lb(ElementType(options.get_double("primal_regularization_lb"))),
         primal_regularization_decrease_factor(ElementType(options.get_double("primal_regularization_decrease_factor"))),
         primal_regularization_fast_increase_factor(ElementType(options.get_double("primal_regularization_fast_increase_factor"))),
         primal_regularization_slow_increase_factor(ElementType(options.get_double("primal_regularization_slow_increase_factor"))),
         threshold_unsuccessful_attempts(options.get_unsigned_int("threshold_unsuccessful_attempts")) {
   }

   template <typename ElementType>
   void PrimalDualRegularization<ElementType>::initialize_memory(const OptimizationProblem& problem, const HessianModel& hessian_model) {
      this->number_variables = problem.number_variables;
      this->number_constraints = problem.number_constraints;
      this->dimension = problem.number_variables + problem.number_constraints;
      this->number_hessian_nonzeros = problem.number_hessian_nonzeros(hessian_model);
      this->number_jacobian_nonzeros = problem.number_jacobian_nonzeros();
   }

   template <typename ElementType>
   void PrimalDualRegularization<ElementType>::initialize_statistics(Statistics& statistics, const Options& options) {
      statistics.add_column("regulariz", Statistics::double_width - 4, options.get_int("statistics_regularization_column_order"));
   }

   template <typename ElementType>
   void PrimalDualRegularization<ElementType>::regularize_hessian(Statistics& statistics, const Subproblem& subproblem,
         SymmetricMatrix<size_t, ElementType>& hessian, const Inertia& expected_inertia) {
      // pick the member linear solver
      if (this->optional_linear_solver == nullptr) {
         this->optional_linear_solver = SymmetricIndefiniteLinearSolverFactory::create(this->optional_linear_solver_name);
         this->optional_linear_solver->initialize_memory(subproblem);
      }
      this->regularize_hessian(statistics, subproblem, hessian, expected_inertia, *this->optional_linear_solver);
   }

   template <typename ElementType>
   void PrimalDualRegularization<ElementType>::regularize_hessian(Statistics& /*statistics*/, const Subproblem& subproblem,
         SymmetricMatrix<size_t, ElementType>& /*hessian*/, const Inertia& /*expected_inertia*/,
         DirectSymmetricIndefiniteLinearSolver<size_t, double>& /*linear_solver*/) {
      // to regularize the Hessian only, call the function for the augmented matrix with no dual part
      // TODO fix
      throw std::runtime_error("PrimalDualRegularization::regularize_hessian not implemented yet");
   }

   // the augmented matrix has been factorized prior to calling this function
   template <typename ElementType>
   void PrimalDualRegularization<ElementType>::regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         SymmetricMatrix<size_t, ElementType>& augmented_matrix, ElementType dual_regularization_parameter,
         const Inertia& expected_inertia) {
      if (this->optional_linear_solver == nullptr) {
         this->optional_linear_solver = SymmetricIndefiniteLinearSolverFactory::create(this->optional_linear_solver_name);
         this->optional_linear_solver->initialize_memory(subproblem);
      }
      this->regularize_augmented_matrix(statistics, subproblem, augmented_matrix, dual_regularization_parameter,
         expected_inertia, *this->optional_linear_solver);
   }

   template <typename ElementType>
   void PrimalDualRegularization<ElementType>::regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         SymmetricMatrix<size_t, ElementType>& augmented_matrix, ElementType dual_regularization_parameter,
         const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<size_t, double>& linear_solver) {
      DEBUG2 << "Original matrix\n" << augmented_matrix << '\n';

      this->primal_regularization = ElementType(0);
      this->dual_regularization = ElementType(0);
      DEBUG << "Testing factorization with regularization factors (0, 0)\n";
      size_t number_attempts = 1;
      DEBUG << "Number of attempts: " << number_attempts << "\n\n";

      // perform the symbolic analysis only once
      if (!this->symbolic_analysis_performed) {
         DEBUG << "Performing symbolic analysis of the indefinite system\n";
         linear_solver.do_symbolic_analysis(augmented_matrix);
         this->symbolic_analysis_performed = true;
      }

      DEBUG << "Performing numerical factorization of the indefinite system\n";
      linear_solver.do_numerical_factorization(augmented_matrix);
      const Inertia estimated_inertia = linear_solver.get_inertia();
      DEBUG << "Expected inertia  " << expected_inertia << '\n';
      DEBUG << "Estimated inertia " << estimated_inertia << '\n';

      if (estimated_inertia == expected_inertia) {
         DEBUG << "The inertia is correct\n";
         statistics.set("regulariz", this->primal_regularization);
         return;
      }

      // set the constraint regularization coefficient
      if (linear_solver.matrix_is_singular()) {
         DEBUG << "Matrix is singular\n";
         this->dual_regularization = this->dual_regularization_fraction * dual_regularization_parameter;
      }
      // set the Hessian regularization coefficient
      if (this->previous_primal_regularization == 0.) {
         this->primal_regularization = this->primal_regularization_initial_factor;
      }
      else {
         this->primal_regularization = std::max(this->primal_regularization_lb,
            this->previous_primal_regularization / this->primal_regularization_decrease_factor);
      }

      // regularize the augmented matrix
      const Collection<size_t>& primal_indices = subproblem.get_primal_regularization_variables();
      const Collection<size_t>& dual_indices = subproblem.get_dual_regularization_constraints();
      // offset 0
      augmented_matrix.set_regularization(primal_indices, 0, this->primal_regularization);
      // offset is the size of the primal indices
      augmented_matrix.set_regularization(dual_indices, primal_indices.size(), -this->dual_regularization);

      bool good_inertia = false;
      while (!good_inertia) {
         DEBUG << "Testing factorization with regularization factors (" << this->primal_regularization << ", " << this->dual_regularization << ")\n";
         DEBUG2 << augmented_matrix << '\n';
         DEBUG << "Performing numerical factorization of the indefinite system\n";
         linear_solver.do_numerical_factorization(augmented_matrix);
         number_attempts++;
         DEBUG << "Number of attempts: " << number_attempts << "\n";

         const Inertia new_estimated_inertia = linear_solver.get_inertia();
         DEBUG << "Expected inertia  " << expected_inertia << '\n';
         DEBUG << "Estimated inertia " << new_estimated_inertia << '\n';

         if (new_estimated_inertia.positive == expected_inertia.positive && new_estimated_inertia.negative == expected_inertia.negative &&
               new_estimated_inertia.zero == expected_inertia.zero) {
            good_inertia = true;
            DEBUG << "The inertia is correct\n";
            this->previous_primal_regularization = this->primal_regularization;
         }
         else {
            if (this->previous_primal_regularization == 0. || this->threshold_unsuccessful_attempts < number_attempts) {
               this->primal_regularization *= this->primal_regularization_fast_increase_factor;
            }
            else {
               this->primal_regularization *= this->primal_regularization_slow_increase_factor;
            }

            if (this->primal_regularization <= this->regularization_failure_threshold) {
               // regularize the augmented matrix
               augmented_matrix.set_regularization(primal_indices, 0, this->primal_regularization);
               augmented_matrix.set_regularization(dual_indices, primal_indices.size(), -this->dual_regularization);
            }
            else {
               throw UnstableRegularization();
            }
         }
      }
      statistics.set("regulariz", this->primal_regularization);
   }

   template <typename ElementType>
   bool PrimalDualRegularization<ElementType>::performs_primal_regularization() const {
      return true;
   }

   template <typename ElementType>
   bool PrimalDualRegularization<ElementType>::performs_dual_regularization() const {
      return true;
   }

   template <typename ElementType>
   [[nodiscard]] double PrimalDualRegularization<ElementType>::get_primal_regularization_factor() const {
      return this->primal_regularization;
   }

   template <typename ElementType>
   std::string PrimalDualRegularization<ElementType>::get_name() const {
      return "primal-dual";
   }
} // namespace

#endif // UNO_PRIMALDUALREGULARIZATION_H