#ifndef CONSTRAINTRELAXATIONSTRATEGY_H
#define CONSTRAINTRELAXATIONSTRATEGY_H

#include <vector>
#include <cmath>
#include "Statistics.hpp"
#include "Subproblem.hpp"
#include "Direction.hpp"
#include "Problem.hpp"
#include "Iterate.hpp"

class ConstraintRelaxationStrategy {
public:
   explicit ConstraintRelaxationStrategy(Subproblem& subproblem);
   virtual Iterate initialize(Statistics& statistics, const Problem& problem, std::vector<double>& x, Multipliers& multipliers) = 0;

   virtual void generate_subproblem(const Problem& problem, const Iterate& current_iterate, double objective_multiplier, double trust_region_radius)
   = 0;
   virtual Direction compute_feasible_direction(const Problem& problem, Iterate& current_iterate) = 0;
   virtual Direction solve_feasibility_problem(const Problem& problem, Iterate& current_iterate, Direction& direction) = 0;
   void update_variables_bounds(const Problem& problem, const Iterate& current_iterate, double trust_region_radius);

   virtual bool is_acceptable(Statistics& statistics, const Problem& problem, Iterate& current_iterate, Iterate& trial_iterate, Direction&
   direction, double step_length) = 0;
   virtual double compute_predicted_reduction(const Problem& problem, Iterate& current_iterate, Direction& direction, double step_length) = 0;

   [[nodiscard]] int get_hessian_evaluation_count() const;
   [[nodiscard]] int get_number_subproblems_solved() const;

protected:
   Subproblem& subproblem;
};

#endif //CONSTRAINTRELAXATIONSTRATEGY_H
