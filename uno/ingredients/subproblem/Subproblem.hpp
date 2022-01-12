#ifndef UNO_SUBPROBLEM_H
#define UNO_SUBPROBLEM_H

#include <vector>
#include <memory>
#include "optimization/Problem.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/Constraint.hpp"
#include "Direction.hpp"
#include "PredictedReductionModel.hpp"
#include "linear_algebra/Vector.hpp"
#include "solvers/linear/LinearSolver.hpp"
#include "tools/Statistics.hpp"

enum SecondOrderCorrection {
   NO_SOC = 0,
   SOC_UPON_ACCEPTANCE,
   SOC_UPON_REJECTION
};

/*! \class Subproblem
 * \brief Subproblem
 *
 *  Local approximation of a nonlinear optimization problem (virtual class) 
 */
class Subproblem {
public:
   Subproblem(size_t number_variables, size_t max_number_variables, size_t number_constraints, bool uses_slacks, SecondOrderCorrection soc_strategy,
         bool is_second_order_method, Norm residual_norm);
   virtual ~Subproblem() = default;

   // virtual methods implemented by subclasses
   virtual void initialize(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& first_iterate);
   virtual void create_current_subproblem(const Problem& problem, const Scaling& scaling, Iterate& current_iterate, double objective_multiplier,
         double trust_region_radius) = 0;
   virtual void build_objective_model(const Problem& problem, const Scaling& scaling, Iterate& current_iterate, double objective_multiplier) = 0;
   virtual void add_elastic_variable(Iterate& current_iterate, size_t i, double objective_term, size_t j, double jacobian_term);
   virtual void remove_elastic_variable(size_t i, size_t j);

   // direction computation
   virtual Direction solve(Statistics& statistics, const Problem& problem, Iterate& current_iterate) = 0;
   virtual Direction compute_second_order_correction(const Problem& problem, Iterate& trial_iterate);

   // globalization metrics
   [[nodiscard]] virtual PredictedReductionModel generate_predicted_reduction_model(const Problem& problem, const Direction& direction) const = 0;
   virtual void compute_progress_measures(const Problem& problem, const Scaling& scaling, Iterate& iterate);
   virtual void register_accepted_iterate(Iterate& iterate);

   [[nodiscard]] virtual size_t get_hessian_evaluation_count() const = 0;
   virtual void set_initial_point(const std::vector<double>& initial_point) = 0;

   // available methods
   void set_scaled_objective_gradient(const Problem& problem, const Scaling& scaling, Iterate& current_iterate, double objective_multiplier);
   // feasibility subproblem
   void compute_feasibility_linear_objective(const Iterate& current_iterate, const ConstraintPartition& constraint_partition);
   void generate_feasibility_bounds(const Problem& problem, const std::vector<double>& current_constraints, const ConstraintPartition&
   constraint_partition);
   static double push_variable_to_interior(double variable_value, const Range& variable_bounds);
   void set_constraints_bounds(const Problem& problem, const std::vector<double>& current_constraints);

   double compute_first_order_error(const Problem& problem, const Scaling& scaling, Iterate& iterate, double objective_multiplier) const;
   void compute_residuals(const Problem& problem, const Scaling& scaling, Iterate& iterate, double objective_multiplier) const;

   static double compute_complementarity_error(const Problem& problem, const Scaling& scaling, Iterate& iterate,
         const std::vector<double>& constraint_multipliers, const std::vector<double>& lower_bounds_multipliers,
         const std::vector<double>& upper_bounds_multipliers);

   size_t number_variables; // can be updated on the fly (elastic variables)
   const size_t max_number_variables;
   const size_t number_constraints;
   const bool uses_slacks;
   const SecondOrderCorrection soc_strategy;
   // when the subproblem is reformulated (e.g. when slacks are introduced), the bounds may be altered
   std::vector<Range> variables_bounds;
   std::vector<double> constraints_multipliers;
   SparseVector<double> objective_gradient;
   std::vector<SparseVector<double>> constraints_jacobian;
   std::vector<Range> constraints_bounds;
   // Hessian is optional and depends on the subproblem
   Direction direction;

   size_t number_subproblems_solved{0};
   // when the parameterization of the subproblem (e.g. penalty or barrier parameter) is updated, signal it
   bool subproblem_definition_changed{false};
   const bool is_second_order_method;
   const Norm residual_norm;

protected:
   virtual void set_variables_bounds(const Problem& problem, const Iterate& current_iterate, double trust_region_radius);
};

#endif // UNO_SUBPROBLEM_H
