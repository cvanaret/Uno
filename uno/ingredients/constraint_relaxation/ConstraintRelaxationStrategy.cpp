#include "ConstraintRelaxationStrategy.hpp"

ConstraintRelaxationStrategy::ConstraintRelaxationStrategy(const Problem& problem, Subproblem& subproblem):
      subproblem(subproblem),
      elastic_variables(ConstraintRelaxationStrategy::count_elastic_variables(problem)),
      // save the original number of variables in the subproblem
      number_subproblem_variables(subproblem.number_variables) {
   // generate elastic variables to relax the constraints
   ConstraintRelaxationStrategy::generate_elastic_variables(problem, this->elastic_variables, subproblem.number_variables);
}

size_t ConstraintRelaxationStrategy::count_elastic_variables(const Problem& problem) {
   size_t number_variables = 0;
   for (size_t j = 0; j < problem.number_constraints; j++) {
      if (-std::numeric_limits<double>::infinity() < problem.constraint_bounds[j].lb) {
         number_variables++;
      }
      if (problem.constraint_bounds[j].ub < std::numeric_limits<double>::infinity()) {
         number_variables++;
      }
   }
   return number_variables;
}

void ConstraintRelaxationStrategy::generate_elastic_variables(const Problem& problem, ElasticVariables& elastic_variables, size_t number_variables) {
   // generate elastic variables p and n on the fly to relax the constraints
   size_t elastic_index = number_variables;
   for (size_t j = 0; j < problem.number_constraints; j++) {
      if (-std::numeric_limits<double>::infinity() < problem.constraint_bounds[j].lb) {
         // nonpositive variable n that captures the negative part of the constraint violation
         elastic_variables.negative.insert(j, elastic_index);
         elastic_index++;
      }
      if (problem.constraint_bounds[j].ub < std::numeric_limits<double>::infinity()) {
         // nonnegative variable p that captures the positive part of the constraint violation
         elastic_variables.positive.insert(j, elastic_index);
         elastic_index++;
      }
   }
}

void ConstraintRelaxationStrategy::add_elastic_variables_to_subproblem() {
   const double objective_coefficient = 1.;
   // add the positive elastic variables
   this->elastic_variables.positive.for_each([&](size_t j, size_t i) {
      this->subproblem.add_variable(i, 0., {0., std::numeric_limits<double>::infinity()}, objective_coefficient, j, -1.);
   });
   this->elastic_variables.negative.for_each([&](size_t j, size_t i) {
      this->subproblem.add_variable(i, 0., {0., std::numeric_limits<double>::infinity()}, objective_coefficient, j, 1.);
   });
}

void ConstraintRelaxationStrategy::remove_elastic_variables_from_subproblem() {
   const auto erase_elastic_variables = [&](size_t j, size_t i) {
      this->subproblem.remove_variable(i, j);
   };
   this->elastic_variables.positive.for_each(erase_elastic_variables);
   this->elastic_variables.negative.for_each(erase_elastic_variables);
}

void ConstraintRelaxationStrategy::remove_elastic_variables_from_direction(const Problem& problem, Direction& direction) {
   // the primal variables and corresponding bound multipliers are organized as follows:
   // original | subproblem-specific (may be empty) | elastic
   direction.x.resize(this->number_subproblem_variables);
   direction.multipliers.lower_bounds.resize(this->number_subproblem_variables);
   direction.multipliers.upper_bounds.resize(this->number_subproblem_variables);
   direction.norm = norm_inf(direction.x);
   // recover active set
   this->recover_active_set(problem, direction);
}

void ConstraintRelaxationStrategy::recover_active_set(const Problem& problem, Direction& direction) {
   // TODO remove elastic variables
   for (size_t i = problem.number_variables; i < direction.x.size(); i++) {
      //direction.active_set.bounds.at_lower_bound.erase(i);
      //direction.active_set.bounds.at_upper_bound.erase(i);
   }
   // constraints: only when p-n = 0
   if (direction.constraint_partition.has_value()) {
      ConstraintPartition& constraint_partition = direction.constraint_partition.value();
      this->elastic_variables.positive.for_each([&](size_t j, size_t i) {
         // if the component is strictly positive, the constraint is violated
         if (0. < direction.x[i]) {
            constraint_partition.lower_bound_infeasible.push_back(j);
         }
      });
      this->elastic_variables.negative.for_each([&](size_t j, size_t i) {
         // if the component is strictly positive, the constraint is violated
         if (0. < direction.x[i]) {
            constraint_partition.upper_bound_infeasible.push_back(j);
         }
      });
   }

   /*
   for (size_t j = 0; j < direction.multipliers.constraints.size(); j++) {
      // compute constraint violation
      double constraint_violation = 0.;
      try {
         size_t i = this->elastic_variables.positive.at(j);
         constraint_violation += direction.x[i];
      }
      catch (const std::out_of_range& e) {
      }
      try {
         size_t i = this->elastic_variables.negative.at(j);
         constraint_violation += direction.x[i];
      }
      catch (const std::out_of_range& e) {
      }
      // update active set
      if (0. < constraint_violation) {
         //direction.active_set.constraints.at_lower_bound.erase(j);
         //direction.active_set.constraints.at_upper_bound.erase(j);
      }
   }
    */
}

Direction ConstraintRelaxationStrategy::compute_second_order_correction(const Problem& problem, Iterate& trial_iterate) {
   return this->subproblem.compute_second_order_correction(problem, trial_iterate);
}

PredictedReductionModel ConstraintRelaxationStrategy::generate_predicted_reduction_model(const Problem& problem, const Direction& direction) const {
   return this->subproblem.generate_predicted_reduction_model(problem, direction);
}

size_t ConstraintRelaxationStrategy::get_hessian_evaluation_count() const {
   return this->subproblem.get_hessian_evaluation_count();
}

size_t ConstraintRelaxationStrategy::get_number_subproblems_solved() const {
   return this->subproblem.number_subproblems_solved;
}

SecondOrderCorrection ConstraintRelaxationStrategy::soc_strategy() const {
   return this->subproblem.soc_strategy;
}

void ConstraintRelaxationStrategy::register_accepted_iterate(Iterate& iterate) {
   this->subproblem.register_accepted_iterate(iterate);
}