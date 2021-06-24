#include "ConstraintRelaxationStrategy.hpp"

ConstraintRelaxationStrategy::ConstraintRelaxationStrategy(std::unique_ptr<Subproblem> subproblem): subproblem(std::move(subproblem)) {
}

void ConstraintRelaxationStrategy::update_variables_bounds(const Problem& problem, const Iterate& current_iterate, double trust_region_radius) {
   this->subproblem->set_trust_region(problem, current_iterate, trust_region_radius);
}

int ConstraintRelaxationStrategy::get_hessian_evaluation_count() const {
   return this->subproblem->get_hessian_evaluation_count();
}

int ConstraintRelaxationStrategy::get_number_subproblems_solved() const {
   return this->subproblem->number_subproblems_solved;
}