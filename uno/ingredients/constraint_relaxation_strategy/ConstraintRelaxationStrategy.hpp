// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CONSTRAINTRELAXATIONSTRATEGY_H
#define UNO_CONSTRAINTRELAXATIONSTRATEGY_H

#include "ingredients/globalization_strategy/ProgressMeasures.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/Statistics.hpp"
#include "tools/Options.hpp"

class ConstraintRelaxationStrategy {
public:
   ConstraintRelaxationStrategy(const Model& model, const Options& options);
   virtual ~ConstraintRelaxationStrategy() = default;

   virtual void initialize(Iterate& initial_iterate) = 0;
   virtual void set_trust_region_radius(double trust_region_radius) = 0;

   // direction computation
   [[nodiscard]] virtual Direction compute_feasible_direction(Statistics& statistics, Iterate& current_iterate,
         WarmstartInformation& warmstart_information) = 0;
   [[nodiscard]] virtual Direction compute_feasible_direction(Statistics& statistics, Iterate& current_iterate,
         const std::vector<double>& initial_point, WarmstartInformation& warmstart_information) = 0;
   virtual void switch_to_feasibility_problem(Iterate& current_iterate, WarmstartInformation& warmstart_information) = 0;

   // trial iterate acceptance
   virtual void compute_progress_measures(Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length) = 0;
   [[nodiscard]] virtual bool is_iterate_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length) = 0;

   [[nodiscard]] virtual size_t get_hessian_evaluation_count() const = 0;
   [[nodiscard]] virtual size_t get_number_subproblems_solved() const = 0;

protected:
   const Model& original_model;
   const Norm progress_norm;
   const Norm residual_norm;
   const double residual_scaling_threshold;

   void compute_primal_dual_residuals(const Model& model, const RelaxedProblem& feasibility_problem, Iterate& iterate);
   static void evaluate_lagrangian_gradient(size_t number_variables, Iterate& iterate, const Multipliers& multipliers, double objective_multiplier);

   [[nodiscard]] double compute_stationarity_error(const Iterate& iterate) const;
   [[nodiscard]] virtual double compute_complementarity_error(const std::vector<double>& primals, const std::vector<double>& constraints,
         const Multipliers& multipliers) const = 0;
   [[nodiscard]] double compute_stationarity_scaling(const Model& model, const Iterate& iterate) const;
   [[nodiscard]] double compute_complementarity_scaling(const Model& model, const Iterate& iterate) const;
};

#endif //UNO_CONSTRAINTRELAXATIONSTRATEGY_H
