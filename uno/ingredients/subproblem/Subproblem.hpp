// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SUBPROBLEM_H
#define UNO_SUBPROBLEM_H

#include <vector>
#include <memory>
#include "optimization/Model.hpp"
#include "optimization/Iterate.hpp"
#include "reformulation/l1RelaxedProblem.hpp"
#include "Direction.hpp"
#include "PredictedOptimalityReductionModel.hpp"
#include "linear_algebra/Vector.hpp"
#include "solvers/linear/LinearSolver.hpp"
#include "tools/Statistics.hpp"

/*! \class Subproblem
 * \brief Subproblem
 */
class Subproblem {
public:
   Subproblem(size_t max_number_variables, size_t max_number_constraints);
   virtual ~Subproblem() = default;

   // virtual methods implemented by subclasses
   virtual void initialize(Statistics& statistics, const NonlinearProblem& problem, Iterate& first_iterate) = 0;
   virtual Direction solve(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate) = 0;
   virtual Direction compute_second_order_correction(const NonlinearProblem& model, Iterate& trial_iterate) = 0;

   void set_variable_bounds(const NonlinearProblem& problem, const Iterate& current_iterate, double trust_region_radius);
   virtual void prepare_for_feasibility_problem(const NonlinearProblem& problem, Iterate& current_iterate) = 0;
   [[nodiscard]] virtual double get_proximal_coefficient() const = 0;
   virtual void set_elastic_variables(const l1RelaxedProblem& problem, Iterate& current_iterate) = 0;

   // globalization metrics
   [[nodiscard]] virtual PredictedOptimalityReductionModel generate_predicted_optimality_reduction_model(const NonlinearProblem& problem,
         const Direction& direction) const = 0;
   [[nodiscard]] virtual double compute_optimality_measure(const NonlinearProblem& problem, Iterate& iterate) = 0;

   virtual void postprocess_accepted_iterate(const NonlinearProblem& model, Iterate& iterate) = 0;

   [[nodiscard]] virtual size_t get_hessian_evaluation_count() const = 0;
   virtual void set_initial_point(const std::vector<double>& initial_point) = 0;

   SparseVector<double> objective_gradient;
   std::vector<double> constraints;
   std::vector<SparseVector<double>> constraint_jacobian;
   std::vector<Interval> variable_bounds;
   Direction direction;

   size_t number_subproblems_solved{0};
   // when the parameterization of the subproblem (e.g. penalty or barrier parameter) is updated, signal it
   bool subproblem_definition_changed{false};

protected:
   static void check_unboundedness(const Direction& direction);
};

#endif // UNO_SUBPROBLEM_H
