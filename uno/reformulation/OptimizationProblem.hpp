// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_OPTIMIZATIONPROBLEM_H
#define UNO_OPTIMIZATIONPROBLEM_H

#include <vector>
#include "optimization/Iterate.hpp"
#include "optimization/Model.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "ingredients/subproblem/Direction.hpp"

class OptimizationProblem {
public:
   OptimizationProblem(const Model& model, size_t number_variables, size_t number_constraints);
   virtual ~OptimizationProblem() = default;

   const Model& model;
   const size_t number_variables; /*!< Number of variables */
   const size_t number_constraints; /*!< Number of constraints */

   [[nodiscard]] bool is_constrained() const;
   [[nodiscard]] bool has_inequality_constraints() const;

   // function evaluations
   [[nodiscard]] virtual double get_objective_multiplier() const = 0;
   virtual void evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const = 0;
   virtual void evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const = 0;
   virtual void evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const = 0;
   virtual void evaluate_lagrangian_hessian(const std::vector<double>& x, const std::vector<double>& multipliers, SymmetricMatrix<double>& hessian) const = 0;

   /*
   virtual void set_infeasibility_measure(Iterate& iterate, Norm progress_norm) const = 0;
   virtual void set_objective_measure(Iterate& iterate) const = 0;
   [[nodiscard]] virtual double compute_predicted_infeasibility_reduction_model(const Iterate& current_iterate, const Direction& direction,
         double step_length, Norm progress_norm) const = 0;
   [[nodiscard]] virtual std::function<double(double)> compute_predicted_objective_reduction_model(const Iterate& current_iterate,
         const Direction& direction, double step_length, const SymmetricMatrix<double>& hessian) const = 0;
*/

   [[nodiscard]] size_t get_number_original_variables() const;
   [[nodiscard]] virtual double variable_lower_bound(size_t variable_index) const = 0;
   [[nodiscard]] virtual double variable_upper_bound(size_t variable_index) const = 0;
   [[nodiscard]] virtual double constraint_lower_bound(size_t constraint_index) const = 0;
   [[nodiscard]] virtual double constraint_upper_bound(size_t constraint_index) const = 0;
   [[nodiscard]] virtual const Collection<size_t>& get_lower_bounded_variables() const = 0;
   [[nodiscard]] virtual const Collection<size_t>& get_upper_bounded_variables() const = 0;
   [[nodiscard]] virtual const Collection<size_t>& get_single_lower_bounded_variables() const = 0;
   [[nodiscard]] virtual const Collection<size_t>& get_single_upper_bounded_variables() const = 0;

   [[nodiscard]] virtual size_t number_objective_gradient_nonzeros() const = 0;
   [[nodiscard]] virtual size_t number_jacobian_nonzeros() const = 0;
   [[nodiscard]] virtual size_t number_hessian_nonzeros() const = 0;
};

inline OptimizationProblem::OptimizationProblem(const Model& model, size_t number_variables, size_t number_constraints):
      model(model), number_variables(number_variables), number_constraints(number_constraints) {
}

inline bool OptimizationProblem::is_constrained() const {
   return (0 < this->number_constraints);
}

inline bool OptimizationProblem::has_inequality_constraints() const {
   return (not this->model.get_inequality_constraints().is_empty());
}

inline size_t OptimizationProblem::get_number_original_variables() const {
   return this->model.number_variables;
}

#endif // UNO_OPTIMIZATIONPROBLEM_H
