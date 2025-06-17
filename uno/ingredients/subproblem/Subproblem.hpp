// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SUBPROBLEM_H
#define UNO_SUBPROBLEM_H
#include "ingredients/constraint_relaxation_strategies/OptimizationProblem.hpp"
#include "ingredients/regularization_strategies/RegularizationStrategy.hpp"
#include "optimization/Multipliers.hpp"
#include "optimization/WarmstartInformation.hpp"

namespace uno {
   // forward declarations
   class HessianModel;
   class Iterate;
   class Multipliers;
   template <typename ElementType>
   class RectangularMatrix;
   template <typename ElementType>
   class RegularizationStrategy;
   template <typename ElementType>
   class SparseVector;

   class Subproblem {
   public:
      const size_t number_variables, number_constraints;

      Subproblem(const OptimizationProblem& problem, Iterate& current_iterate, const Multipliers& current_multipliers,
         HessianModel& hessian_model, RegularizationStrategy<double>& regularization_strategy, double trust_region_radius);

      void evaluate_functions(Statistics& statistics, SparseVector<double>& linear_objective,
         std::vector<double>& constraints, RectangularMatrix<double>& constraint_jacobian,
         SymmetricMatrix<size_t, double>& hessian, const WarmstartInformation& warmstart_information) const;

      void regularize_hessian(Statistics& statistics, SymmetricMatrix<size_t, double>& hessian,
         const WarmstartInformation& warmstart_information) const;

      void set_variables_bounds(std::vector<double>& variables_lower_bounds, std::vector<double>& variables_upper_bounds,
         const WarmstartInformation& warmstart_information) const;

      template <typename Array>
      void set_constraints_bounds(Array& constraints_lower_bounds, Array& constraints_upper_bounds,
         std::vector<double>& constraints, const WarmstartInformation& warmstart_information) const;

   protected:
      const OptimizationProblem& problem;
      Iterate& current_iterate;
      const Multipliers& current_multipliers;
      HessianModel& hessian_model;
      RegularizationStrategy<double>& regularization_strategy;
      const double trust_region_radius;
   };

   template <typename Array>
   void Subproblem::set_constraints_bounds(Array& constraints_lower_bounds, Array& constraints_upper_bounds,
         std::vector<double>& constraints, const WarmstartInformation& warmstart_information) const {
      if (warmstart_information.constraint_bounds_changed || warmstart_information.constraints_changed) {
         for (size_t constraint_index: Range(this->problem.number_constraints)) {
            constraints_lower_bounds[constraint_index] = this->problem.constraint_lower_bound(constraint_index) - constraints[constraint_index];
            constraints_upper_bounds[constraint_index] = this->problem.constraint_upper_bound(constraint_index) - constraints[constraint_index];
         }
      }
   }
} // namespace

#endif // UNO_SUBPROBLEM_H