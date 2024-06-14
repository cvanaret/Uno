// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CONSTRAINTRELAXATIONSTRATEGY_H
#define UNO_CONSTRAINTRELAXATIONSTRATEGY_H

#include <cstddef>
#include <memory>
#include "ingredients/globalization_strategy/GlobalizationStrategy.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/Norm.hpp"

// forward declarations
class Direction;
class Iterate;
class Model;
struct Multipliers;
class OptimizationProblem;
class Options;
class Statistics;
template <typename ElementType>
class SymmetricMatrix;
template <typename ElementType>
class Vector;
struct WarmstartInformation;

class ConstraintRelaxationStrategy {
public:
   ConstraintRelaxationStrategy(const Model& model, size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros,
         size_t number_jacobian_nonzeros, size_t number_hessian_nonzeros, const Options& options);
   virtual ~ConstraintRelaxationStrategy() = default;

   virtual void initialize(Statistics& statistics, Iterate& initial_iterate, const Options& options) = 0;
   void set_trust_region_radius(double trust_region_radius);

   [[nodiscard]] virtual size_t maximum_number_variables() const = 0;
   [[nodiscard]] virtual size_t maximum_number_constraints() const = 0;

   // direction computation
   virtual void compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, Direction& direction,
         WarmstartInformation& warmstart_information) = 0;
   void compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, Direction& direction,
         const Vector<double>& initial_point, WarmstartInformation& warmstart_information);
   [[nodiscard]] virtual bool solving_feasibility_problem() const = 0;
   virtual void switch_to_feasibility_problem(Statistics& statistics, Iterate& current_iterate) = 0;

   // trial iterate acceptance
   [[nodiscard]] virtual bool is_iterate_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length) = 0;

   // primal-dual residuals
   virtual void compute_primal_dual_residuals(Iterate& iterate) = 0;

   [[nodiscard]] size_t get_hessian_evaluation_count() const;
   [[nodiscard]] size_t get_number_subproblems_solved() const;

protected:
   const Model& model;
   const std::unique_ptr<GlobalizationStrategy> globalization_strategy;
   const std::unique_ptr<Subproblem> subproblem;
   const Norm progress_norm;
   const Norm residual_norm;
   const double residual_scaling_threshold;

   void set_objective_measure(Iterate& iterate) const;
   void set_infeasibility_measure(Iterate& iterate) const;
   [[nodiscard]] double compute_predicted_infeasibility_reduction_model(const Iterate& current_iterate, const Vector<double>& primal_direction,
         double step_length) const;
   [[nodiscard]] std::function<double(double)> compute_predicted_objective_reduction_model(const Iterate& current_iterate,
         const Vector<double>& primal_direction, double step_length, const SymmetricMatrix<double>& hessian) const;
   void compute_progress_measures(Iterate& current_iterate, Iterate& trial_iterate);
   virtual void evaluate_progress_measures(Iterate& iterate) const = 0;

   void compute_primal_dual_residuals(const OptimizationProblem& optimality_problem, const OptimizationProblem& feasibility_problem, Iterate& iterate);
   void evaluate_lagrangian_gradient(Iterate& iterate, const Multipliers& multipliers) const;

   [[nodiscard]] double compute_stationarity_scaling(const Multipliers& multipliers) const;
   [[nodiscard]] double compute_complementarity_scaling(const Multipliers& multipliers) const;

   void set_statistics(Statistics& statistics, const Iterate& iterate) const;
   void set_progress_statistics(Statistics& statistics, const Iterate& iterate) const;
   virtual void set_dual_residuals_statistics(Statistics& statistics, const Iterate& iterate) const = 0;
};

#endif //UNO_CONSTRAINTRELAXATIONSTRATEGY_H
