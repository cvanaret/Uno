// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NONLINEARPROBLEM_H
#define UNO_NONLINEARPROBLEM_H

#include <string>
#include <vector>
#include <map>
#include "optimization/Iterate.hpp"
#include "optimization/Model.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"
#include "ingredients/subproblem/Direction.hpp"

enum Smoothness {
   SMOOTH = 0,
   NONSMOOTH
};

class ReformulatedProblem {
public:
   ReformulatedProblem(const Model& model, size_t number_variables, size_t number_constraints);
   virtual ~ReformulatedProblem() = default;

   const Model& model;
   const size_t number_variables; /*!< Number of variables */
   const size_t number_constraints; /*!< Number of constraints */

   [[nodiscard]] bool is_constrained() const;

   SparseVector<size_t> equality_constraints; /*!< inequality constraints */
   SparseVector<size_t> inequality_constraints; /*!< inequality constraints */
   // lists of bounded variables
   std::vector<size_t> lower_bounded_variables{}; // indices of the lower-bounded variables
   std::vector<size_t> upper_bounded_variables{}; // indices of the upper-bounded variables

   // function evaluations
   [[nodiscard]] virtual double get_objective_multiplier() const = 0;
   [[nodiscard]] virtual double evaluate_objective(Iterate& iterate) const = 0;
   virtual void evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const = 0;
   virtual void evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const = 0;
   virtual void evaluate_constraint_jacobian(Iterate& iterate, std::vector<SparseVector<double>>& constraint_jacobian) const = 0;
   virtual void evaluate_lagrangian_hessian(const std::vector<double>& x, const std::vector<double>& multipliers, SymmetricMatrix& hessian) const = 0;
   //void evaluate_lagrangian_gradient(Iterate& iterate, std::vector<double>& lagrangian_gradient) const;

   [[nodiscard]] virtual double predicted_reduction_contribution(const Iterate& current_iterate, const Direction& direction, double step_length) const = 0;

   [[nodiscard]] size_t get_number_original_variables() const;
   [[nodiscard]] virtual double get_variable_lower_bound(size_t i) const = 0;
   [[nodiscard]] virtual double get_variable_upper_bound(size_t i) const = 0;
   [[nodiscard]] virtual double get_constraint_lower_bound(size_t j) const = 0;
   [[nodiscard]] virtual double get_constraint_upper_bound(size_t j) const = 0;

   [[nodiscard]] virtual size_t get_maximum_number_hessian_nonzeros() const = 0;
};

inline ReformulatedProblem::ReformulatedProblem(const Model& model, size_t number_variables, size_t number_constraints):
      model(model), number_variables(number_variables), number_constraints(number_constraints) {
}

inline bool ReformulatedProblem::is_constrained() const {
   return (0 < this->number_constraints);
}

/*
inline void NonlinearReformulation::evaluate_lagrangian_gradient(Iterate& iterate, std::vector<double>& lagrangian_gradient) const {
   iterate.evaluate_lagrangian_gradient(this->model, iterate.multipliers.constraints, iterate.multipliers.lower_bounds,
         iterate.multipliers.upper_bounds);
   lagrangian_gradient = iterate.lagrangian_gradient;
}
*/

inline size_t ReformulatedProblem::get_number_original_variables() const {
   return this->model.number_variables;
}

#endif // UNO_NONLINEARPROBLEM_H
