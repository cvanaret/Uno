#ifndef CONSTRAINTRELAXATIONSTRATEGY_H
#define CONSTRAINTRELAXATIONSTRATEGY_H

#include <vector>
#include <cmath>
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem/Direction.hpp"
#include "optimization_problem/Problem.hpp"
#include "optimization_problem/Iterate.hpp"
#include "tools/Statistics.hpp"

struct ElasticVariables {
   std::map<size_t, size_t> positive;
   std::map<size_t, size_t> negative;
   [[nodiscard]] size_t size() const { return this->positive.size() + this->negative.size(); }
};

class ConstraintRelaxationStrategy {
public:
   explicit ConstraintRelaxationStrategy(std::unique_ptr<Subproblem> subproblem);
   virtual Iterate initialize(Statistics& statistics, const Problem& problem, std::vector<double>& x, Multipliers& multipliers) = 0;

   virtual void generate_subproblem(const Problem& problem, Iterate& current_iterate, double trust_region_radius) = 0;

   virtual Direction compute_feasible_direction(Statistics& statistics, const Problem& problem, Iterate& current_iterate) = 0;
   virtual Direction solve_feasibility_problem(Statistics& statistics, const Problem& problem, Iterate& current_iterate, Direction& direction) = 0;

   Direction compute_second_order_correction(const Problem& problem, Iterate& trial_iterate);

   virtual bool is_acceptable(Statistics& statistics, const Problem& problem, Iterate& current_iterate, Iterate& trial_iterate, const Direction&
   direction, double step_length) = 0;
   virtual double compute_predicted_reduction(const Problem& problem, Iterate& current_iterate, const Direction& direction, double step_length) = 0;
   virtual void register_accepted_iterate(Iterate& iterate);

   [[nodiscard]] int get_hessian_evaluation_count() const;
   [[nodiscard]] int get_number_subproblems_solved() const;

protected:
   std::unique_ptr<Subproblem> subproblem;

   static void generate_elastic_variables(const Problem& problem, ElasticVariables& elastic_variables);
   void set_elastic_bounds_in_subproblem(const Problem& problem, size_t number_elastic_variables) const;
   void add_elastic_variables_to_subproblem(const ElasticVariables& elastic_variables);
};

#endif //CONSTRAINTRELAXATIONSTRATEGY_H
