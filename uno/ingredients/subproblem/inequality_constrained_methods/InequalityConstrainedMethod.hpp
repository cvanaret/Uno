// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INEQUALITYCONSTRAINEDMETHOD_H
#define UNO_INEQUALITYCONSTRAINEDMETHOD_H

#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   class InequalityConstrainedMethod : public Subproblem {
   public:
      InequalityConstrainedMethod(size_t number_variables, size_t number_constraints);
      ~InequalityConstrainedMethod() override = default;
      
      void initialize_statistics(Statistics& statistics, const Options& options) override;
      void set_initial_point(const Vector<double>& point) override;
      void initialize_feasibility_problem(const l1RelaxedProblem& problem, Iterate& current_iterate) override;
      void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate) override;
      void exit_feasibility_problem(const OptimizationProblem& problem, Iterate& trial_iterate) override;

      void set_auxiliary_measure(const Model& model, Iterate& iterate) override;
      [[nodiscard]] double compute_predicted_auxiliary_reduction_model(const Model& model, const Iterate&, const Vector<double>&, double) const override;

      void postprocess_iterate(const OptimizationProblem& model, Iterate& iterate) override;

   protected:
      Vector<double> initial_point{};
      std::vector<double> direction_lower_bounds{};
      std::vector<double> direction_upper_bounds{};
      std::vector<double> linearized_constraints_lower_bounds{};
      std::vector<double> linearized_constraints_upper_bounds{};

      SparseVector<double> objective_gradient; /*!< Sparse Jacobian of the objective */
      std::vector<double> constraints; /*!< Constraint values (size \f$m)\f$ */
      RectangularMatrix<double> constraint_jacobian; /*!< Sparse Jacobian of the constraints */

      void set_direction_bounds(const OptimizationProblem& problem, const Iterate& current_iterate);
      void set_linearized_constraint_bounds(const OptimizationProblem& problem, const std::vector<double>& current_constraints);
      static void compute_dual_displacements(const Multipliers& current_multipliers, Multipliers& direction_multipliers);
   };
} // namespace

#endif // UNO_INEQUALITYCONSTRAINEDMETHOD_H
