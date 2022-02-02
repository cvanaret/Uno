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
 */
class Subproblem {
public:
   Subproblem(size_t max_number_variables, size_t number_constraints, SecondOrderCorrection soc_strategy, bool is_second_order_method,
         Norm residual_norm);
   virtual ~Subproblem() = default;

   // virtual methods implemented by subclasses
   virtual void initialize(Statistics& statistics, const Problem& problem, Iterate& first_iterate);
   virtual void build_current_subproblem(const Problem& problem, Iterate& current_iterate, double objective_multiplier,
         double trust_region_radius) = 0;
   virtual void build_objective_model(const Problem& problem, Iterate& current_iterate, double objective_multiplier) = 0;
   [[nodiscard]] virtual double get_proximal_coefficient() const = 0;

   // direction computation
   virtual Direction solve(Statistics& statistics, const Problem& problem, Iterate& current_iterate) = 0;
   virtual Direction compute_second_order_correction(const Problem& problem, Iterate& trial_iterate);

   // globalization metrics
   [[nodiscard]] virtual PredictedReductionModel generate_predicted_reduction_model(const Problem& problem, const Iterate& current_iterate,
         const Direction& direction) const = 0;
   [[nodiscard]] double compute_first_order_error(const Problem& problem, Iterate& iterate) const;
   [[nodiscard]] virtual double compute_optimality_measure(const Problem& problem, Iterate& iterate);
   [[nodiscard]] static double compute_complementarity_error(const Problem& problem, const Iterate& iterate, const std::vector<double>& constraint_multipliers,
         const std::vector<double>& lower_bounds_multipliers, const std::vector<double>& upper_bounds_multipliers);
   void compute_nonlinear_residuals(const Problem& problem, Iterate& iterate) const;

   virtual void register_accepted_iterate(const Problem& problem, Iterate& iterate);

   [[nodiscard]] virtual size_t get_hessian_evaluation_count() const = 0;
   virtual void set_initial_point(const std::vector<double>& initial_point) = 0;

   void set_scaled_objective_gradient(const Problem& problem, Iterate& current_iterate, double objective_multiplier);
   [[nodiscard]] static double push_variable_to_interior(double variable_value, const Range& variable_bounds);
   void set_constraint_bounds(const Problem& problem, const std::vector<double>& current_constraints);

   const SecondOrderCorrection soc_strategy;
   std::vector<Range> current_variable_bounds;
   SparseVector<double> objective_gradient;
   std::vector<Range> constraint_bounds;
   Direction direction;

   size_t number_subproblems_solved{0};
   // when the parameterization of the subproblem (e.g. penalty or barrier parameter) is updated, signal it
   bool subproblem_definition_changed{false};
   const bool is_second_order_method;
   const Norm residual_norm;

protected:
   virtual void set_current_variable_bounds(const Problem& problem, const Iterate& current_iterate, double trust_region_radius);
};

#endif // UNO_SUBPROBLEM_H
