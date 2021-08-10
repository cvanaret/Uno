#ifndef SUBPROBLEM_H
#define SUBPROBLEM_H

#include <cmath>
#include <vector>
#include <memory>
#include "Statistics.hpp"
#include "Problem.hpp"
#include "Iterate.hpp"
#include "Direction.hpp"
#include "Constraint.hpp"
#include "MA57Solver.hpp"
#include "Vector.hpp"

/*! \class Subproblem
 * \brief Subproblem
 *
 *  Local approximation of a nonlinear optimization problem (virtual class) 
 */
class Subproblem {
public:
   Subproblem(size_t number_variables, size_t number_constraints);
   virtual ~Subproblem() = default;

   virtual void evaluate_constraints(const Problem& problem, Iterate& iterate) const;
   virtual Iterate generate_initial_iterate(Statistics& statistics, const Problem& problem, std::vector<double>& x, Multipliers& multipliers);
   virtual void generate(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) = 0;
   virtual void set_trust_region(const Problem& problem, const Iterate& current_iterate, double trust_region_radius);
   virtual void update_objective_multiplier(const Problem& problem, const Iterate& current_iterate, double objective_multiplier) = 0;

   // direction computation
   virtual Direction compute_direction(Statistics& statistics, const Problem& problem, Iterate& current_iterate) = 0;
   virtual Direction compute_second_order_correction(const Problem& problem, Iterate& trial_iterate);

   // feasibility subproblem
   void compute_feasibility_linear_objective(const Iterate& current_iterate, const ConstraintPartition& constraint_partition);
   void generate_feasibility_bounds(const Problem& problem, const std::vector<double>& current_constraints, const ConstraintPartition&
   constraint_partition);

   // globalization metrics
   virtual double compute_predicted_reduction(const Direction& direction, double step_length) const = 0;
   virtual void compute_progress_measures(const Problem& problem, Iterate& iterate);
   virtual void register_accepted_iterate(Iterate& iterate);

   static double push_variable_to_interior(double variable_value, const Range& variable_bounds);
   void set_constraints_bounds(const Problem& problem, const std::vector<double>& current_constraints);

   static void compute_least_square_multipliers(const Problem& problem, Iterate& current_iterate, std::vector<double>& multipliers, LinearSolver& solver,
         double multipliers_max_size = 1e3);
   static void compute_least_square_multipliers(const Problem& problem, Iterate& current_iterate, std::vector<double>& multipliers,
         double multipliers_max_size = 1e3);

   virtual double compute_constraint_violation(const Problem& problem, const Iterate& iterate) const;
   static double compute_first_order_error(const Problem& problem, Iterate& iterate, double objective_multiplier);
   void compute_errors(const Problem& problem, Iterate& iterate, double objective_multiplier) const;
   virtual int get_hessian_evaluation_count() const = 0;
   virtual void set_initial_point(const std::vector<double>& initial_point) = 0;
   static double compute_complementarity_error(const Problem& problem, Iterate& iterate, const Multipliers& multipliers);

   const size_t number_variables;
   const size_t number_constraints;
   // when the subproblem is reformulated (e.g. when slacks are introduced), the bounds may be altered
   std::vector<Range> variables_bounds;
   std::vector<double> constraints_multipliers;
   SparseVector objective_gradient;
   std::vector<SparseVector> constraints_jacobian;
   std::vector<Range> constraints_bounds;
   // Hessian is optional and depends on the subproblem

   int number_subproblems_solved;
   // when the parameterization of the subproblem (e.g. penalty or barrier parameter) is updated, signal it
   bool subproblem_definition_changed;
   //bool scale_residuals;

protected:
   virtual void set_variables_bounds(const Problem& problem, const Iterate& current_iterate, double trust_region_radius);
};

#endif // SUBPROBLEM_H
