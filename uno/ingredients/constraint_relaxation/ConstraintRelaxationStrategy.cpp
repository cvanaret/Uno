#include "ConstraintRelaxationStrategy.hpp"

ConstraintRelaxationStrategy::ConstraintRelaxationStrategy(const Problem& problem, Subproblem& subproblem):
      subproblem(subproblem),
      //number_variables(this->subproblem.number_variables),
      //number_constraints(this->subproblem.number_constraints),
      elastic_variables(ConstraintRelaxationStrategy::count_elastic_variables(problem)) {
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

void ConstraintRelaxationStrategy::add_elastic_variables_to_subproblem(const ElasticVariables& elastic_variables) {
   // add the positive elastic variables
   elastic_variables.positive.for_each([&](size_t j, size_t i) {
      this->subproblem.add_variable(i, 0., {0., std::numeric_limits<double>::infinity()}, 1., j, -1.);
   });
   elastic_variables.negative.for_each([&](size_t j, size_t i) {
      this->subproblem.add_variable(i, 0., {0., std::numeric_limits<double>::infinity()}, 1., j, 1.);
   });
}

Direction ConstraintRelaxationStrategy::compute_second_order_correction(const Problem& problem, Iterate& trial_iterate) {
   return this->subproblem.compute_second_order_correction(problem, trial_iterate);
}

PredictedReductionModel ConstraintRelaxationStrategy::generate_predicted_reduction_model(const Problem& problem, const Direction& direction) const {
   return this->subproblem.generate_predicted_reduction_model(problem, direction);
}

int ConstraintRelaxationStrategy::get_hessian_evaluation_count() const {
   return this->subproblem.get_hessian_evaluation_count();
}

int ConstraintRelaxationStrategy::get_number_subproblems_solved() const {
   return this->subproblem.number_subproblems_solved;
}

SecondOrderCorrection ConstraintRelaxationStrategy::soc_strategy() const {
   return this->subproblem.soc_strategy;
}

void ConstraintRelaxationStrategy::register_accepted_iterate(Iterate& iterate) {
   this->subproblem.register_accepted_iterate(iterate);
}