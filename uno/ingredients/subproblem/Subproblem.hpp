#ifndef UNO_SUBPROBLEM_H
#define UNO_SUBPROBLEM_H

#include <vector>
#include <memory>
#include "optimization/Problem.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/Constraint.hpp"
#include "optimization/l1ElasticReformulation.hpp"
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
   virtual void evaluate_objective_gradient(const Problem& problem, Iterate& current_iterate);
   virtual void evaluate_constraint_jacobian(const Problem& problem, Iterate& current_iterate);

   virtual void build_objective_model(const Problem& problem, Iterate& current_iterate, double objective_multiplier) = 0;
   virtual void build_constraint_model(const Problem& problem, Iterate& current_iterate) = 0;

   void set_variable_bounds(const Problem& problem, const Iterate& current_iterate, double trust_region_radius);
   [[nodiscard]] virtual double get_proximal_coefficient() const = 0;
   virtual void set_elastic_variables(const l1ElasticReformulation& problem, Iterate& current_iterate) = 0;

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

   virtual void postprocess_accepted_iterate(const Problem& problem, Iterate& iterate);

   [[nodiscard]] virtual size_t get_hessian_evaluation_count() const = 0;
   virtual void set_initial_point(const std::optional<std::vector<double>>& optional_initial_point) = 0;

   const SecondOrderCorrection soc_strategy;
   std::vector<Range> variable_bounds;
   Direction direction;

   size_t number_subproblems_solved{0};
   // when the parameterization of the subproblem (e.g. penalty or barrier parameter) is updated, signal it
   bool subproblem_definition_changed{false};
   const bool is_second_order_method;
   const Norm residual_norm;
};

#endif // UNO_SUBPROBLEM_H
