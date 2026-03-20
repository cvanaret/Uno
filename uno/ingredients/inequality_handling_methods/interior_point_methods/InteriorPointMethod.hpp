// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INTERIORPOINTMETHOD_H
#define UNO_INTERIORPOINTMETHOD_H

#include <memory>
#include "../InequalityHandlingMethod.hpp"
#include "BarrierParameterUpdateStrategy.hpp"
#include "InteriorPointParameters.hpp"
#include "ingredients/constraint_relaxation_strategies/relaxed_problems/l1RelaxedProblem.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/EvaluationCache.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/Parameterization.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   template <typename BarrierProblem>
   class InteriorPointMethod : public InequalityHandlingMethod {
   public:
      explicit InteriorPointMethod(const Options& options);

      void check_problem(const OptimizationProblem& problem, bool uses_trust_region) override;
      void initialize_statistics(Statistics& statistics) override;
      [[nodiscard]] std::unique_ptr<OptimizationProblem> reformulate(const OptimizationProblem& problem,
         Parameterization& parameterization) override;
      [[nodiscard]] bool update_parameterization(Statistics& statistics, const OptimizationProblem& problem,
         const Iterate& current_iterate, Parameterization& parameterization) override;

      void initialize_feasibility_problem(Iterate& current_iterate) override;
      void set_elastic_variable_values(const l1RelaxedProblem& feasibility_problem, Iterate& iterate, Evaluations& evaluations) override;
      [[nodiscard]] double proximal_coefficient() const override;

      [[nodiscard]] std::string get_name() const override;

   protected:
      BarrierParameterUpdateStrategy<BarrierProblem> barrier_parameter_update_strategy;
      double previous_barrier_parameter;
      const InteriorPointParameters parameters;
      const double least_square_multiplier_max_norm;
      const double l1_constraint_violation_coefficient; // (rho in Section 3.3.1 in IPOPT paper)

      bool first_feasibility_iteration{false};

      [[nodiscard]] double barrier_parameter() const;
      //[[nodiscard]] bool is_small_step(const Vector<double>& current_primals, const Vector<double>& direction_primals) const;
   };

   // class template implementation

   template <typename BarrierProblem>
   InteriorPointMethod<BarrierProblem>::InteriorPointMethod(const Options& options):
         InequalityHandlingMethod(),
         barrier_parameter_update_strategy(options),
         previous_barrier_parameter(options.get_double("barrier_initial_parameter")),
         parameters({
               options.get_double("barrier_tau_min"),
               options.get_double("barrier_k_sigma"),
               options.get_double("barrier_regularization_exponent"),
               options.get_double("barrier_small_direction_factor"),
               options.get_double("barrier_push_variable_to_interior_k1"),
               options.get_double("barrier_push_variable_to_interior_k2"),
               options.get_double("barrier_damping_factor"),
               options.get_double("barrier_default_multiplier")
         }),
         least_square_multiplier_max_norm(options.get_double("least_square_multiplier_max_norm")),
         l1_constraint_violation_coefficient(options.get_double("l1_constraint_violation_coefficient")) {
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::check_problem(const OptimizationProblem& problem, bool uses_trust_region) {
      if (uses_trust_region) {
         throw std::runtime_error("A trust-region radius is not supported yet.");
      }
      if (!problem.get_inequality_constraints().empty()) {
         throw std::runtime_error("The problem has inequality constraints. Create an instance of HomogeneousEqualityConstrainedModel");
      }
      if (!problem.get_fixed_variables().empty()) {
         throw std::runtime_error("The problem has fixed variables. Move them to the set of general constraints.");
      }
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::initialize_statistics(Statistics& statistics) {
      statistics.add_column("Barrier", Statistics::double_width, 2, Statistics::column_order.at("Barrier"));
   }

   template <typename BarrierProblem>
   std::unique_ptr<OptimizationProblem> InteriorPointMethod<BarrierProblem>::reformulate(const OptimizationProblem& problem,
         Parameterization& parameterization) {
      parameterization.set("barrier_parameter", this->barrier_parameter());
      return std::make_unique<BarrierProblem>(problem, this->parameters, parameterization);
   }

   template <typename BarrierProblem>
   bool InteriorPointMethod<BarrierProblem>::update_parameterization(Statistics& statistics, const OptimizationProblem& problem,
         const Iterate& current_iterate, Parameterization& parameterization) {
      bool update = false;
      // possibly update the barrier parameter
      if (!this->first_feasibility_iteration) {
         update = this->barrier_parameter_update_strategy.update_barrier_parameter(problem, current_iterate, current_iterate.residuals);
      }
      else {
         this->first_feasibility_iteration = false;
      }
      parameterization.set("barrier_parameter", this->barrier_parameter());
      statistics.set("Barrier", this->barrier_parameter());
      return update;
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::initialize_feasibility_problem(Iterate& current_iterate) {
      this->first_feasibility_iteration = true;

      // temporarily update the objective multiplier
      this->previous_barrier_parameter = this->barrier_parameter();
      const double new_barrier_parameter = std::max(this->barrier_parameter(), current_iterate.primal_feasibility);
      this->barrier_parameter_update_strategy.set_barrier_parameter(new_barrier_parameter);
      DEBUG << "Barrier parameter mu temporarily updated to " << this->barrier_parameter() << '\n';
   }

   template <typename BarrierProblem>
   // set the elastic variables of the current iterate
   void InteriorPointMethod<BarrierProblem>::set_elastic_variable_values(const l1RelaxedProblem& feasibility_problem,
         Iterate& iterate, Evaluations& evaluations) {
      DEBUG << "IPM: setting the elastic variables and their duals\n";

      for (size_t variable_index: Range(feasibility_problem.number_variables)) {
         const double lower_bound = feasibility_problem.variable_lower_bound(variable_index);
         const double upper_bound = feasibility_problem.variable_upper_bound(variable_index);
         if (is_finite(lower_bound)) {
            iterate.multipliers.lower_bounds[variable_index] = this->parameters.default_multiplier;
         }
         if (is_finite(upper_bound)) {
            iterate.multipliers.upper_bounds[variable_index] = -this->parameters.default_multiplier;
         }
      }

      // c(x) - p + n = 0
      // analytical expression for p and n:
      // (mu_over_rho - jacobian_coefficient*this->barrier_constraints[j] + std::sqrt(radical))/2.
      // where jacobian_coefficient = -1 for p, +1 for n
      // Note: IPOPT uses a '+' sign because they define the Lagrangian as f(x) + \lambda^T c(x)
      evaluations.evaluate_constraints(feasibility_problem.model, iterate.primals);
      const double mu = this->barrier_parameter();
      const auto elastic_setting_function = [&](Iterate& iterate, size_t constraint_index, size_t elastic_index, double jacobian_coefficient) {
         // precomputations
         const double constraint_j = evaluations.constraints[constraint_index];
         const double rho = this->l1_constraint_violation_coefficient;
         const double mu_over_rho = mu / rho;
         const double radical = std::pow(constraint_j, 2) + std::pow(mu_over_rho, 2);
         const double sqrt_radical = std::sqrt(radical);

         iterate.primals[elastic_index] = (mu_over_rho - jacobian_coefficient * constraint_j + sqrt_radical) / 2.;
         iterate.multipliers.lower_bounds[elastic_index] = mu / iterate.primals[elastic_index];
         iterate.multipliers.upper_bounds[elastic_index] = 0.;
         assert(0. < iterate.primals[elastic_index] && "The elastic variable is not strictly positive.");
         assert(0. < iterate.multipliers.lower_bounds[elastic_index] && "The elastic dual is not strictly positive.");
      };
      feasibility_problem.set_elastic_variable_values(iterate, elastic_setting_function);
   }

   template <typename BarrierProblem>
   double InteriorPointMethod<BarrierProblem>::proximal_coefficient() const {
      return std::sqrt(this->barrier_parameter());
   }

   // protected member functions

   template <typename BarrierProblem>
   double InteriorPointMethod<BarrierProblem>::barrier_parameter() const {
      return this->barrier_parameter_update_strategy.get_barrier_parameter();
   }

   /*
   // determine if the direction is a "small direction" (Section 3.9 in IPOPT paper)
   template <typename BarrierProblem>
   bool InteriorPointMethod<BarrierProblem>::is_small_step(const Vector<double>& current_primals, const Vector<double>& direction_primals) const {
      const Range variables_range = Range(this->problem->number_variables);
      const VectorExpression relative_direction_size{variables_range, [&](size_t variable_index) {
         return direction_primals[variable_index] / (1 + std::abs(current_primals[variable_index]));
      }};
      static double machine_epsilon = std::numeric_limits<double>::epsilon();
      return (norm_inf(relative_direction_size) <= this->parameters.small_direction_factor * machine_epsilon);
      return false;
   }
   */

   template <typename BarrierProblem>
   std::string InteriorPointMethod<BarrierProblem>::get_name() const {
      return "primal-dual interior-point method";
   }
} // namespace

#endif // UNO_INTERIORPOINTMETHOD_H