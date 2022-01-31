#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"
#include "optimization/Constraint.hpp"

ConstraintRelaxationStrategy::ConstraintRelaxationStrategy(const Problem& problem, double objective_multiplier, const Options& options):
      // save the original problem
      original_problem(problem),
      // create the relaxed problem by introducing elastic variables
      relaxed_problem(problem, objective_multiplier, stod(options.at("elastic_objective_coefficient")), (options.at("use_proximal_term") == "yes")),
      subproblem(SubproblemFactory::create(this->relaxed_problem, this->relaxed_problem.number_variables, options)),
      number_constraints(this->relaxed_problem.number_constraints) {
}

bool ConstraintRelaxationStrategy::is_small_step(const Direction& direction) {
   // return (direction.norm == 0.);
   const double tolerance = 1e-8;
   const double small_step_factor = 100.;
   return (direction.norm <= tolerance / small_step_factor);
}

void ConstraintRelaxationStrategy::recover_active_set(const Problem& problem, Direction& direction) {
   // TODO remove elastic variables
   for (size_t i = problem.number_variables; i < direction.x.size(); i++) {
      //direction.active_set.bounds.at_lower_bound.erase(i);
      //direction.active_set.bounds.at_upper_bound.erase(i);
   }
   // constraints: only when p-n = 0
   if (direction.constraint_partition.has_value()) {
      /* TODO with this->relaxed_problem
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
       */
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

size_t ConstraintRelaxationStrategy::get_hessian_evaluation_count() const {
   return this->subproblem->get_hessian_evaluation_count();
}

size_t ConstraintRelaxationStrategy::get_number_subproblems_solved() const {
   return this->subproblem->number_subproblems_solved;
}

SecondOrderCorrection ConstraintRelaxationStrategy::soc_strategy() const {
   return this->subproblem->soc_strategy;
}