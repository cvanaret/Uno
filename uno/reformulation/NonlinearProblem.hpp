// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NONLINEARPROBLEM_H
#define UNO_NONLINEARPROBLEM_H

// #include <iomanip>
// std::setprecision(20) <<
#include <string>
#include <vector>
#include <map>
#include "optimization/Iterate.hpp"
#include "optimization/Model.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "ingredients/subproblem/Direction.hpp"
#include "tools/Range.hpp"

class NonlinearProblem {
public:
   NonlinearProblem(const Model& model, size_t number_variables, size_t number_constraints);
   virtual ~NonlinearProblem() = default;

   const Model& model;
   const size_t number_variables; /*!< Number of variables */
   const size_t number_constraints; /*!< Number of constraints */

   [[nodiscard]] bool is_constrained() const;

   SparseVector<size_t> equality_constraints{}; /*!< inequality constraints */
   SparseVector<size_t> inequality_constraints{}; /*!< inequality constraints */
   // lists of bounded variables
   std::vector<size_t> lower_bounded_variables{}; // indices of the lower-bounded variables
   std::vector<size_t> upper_bounded_variables{}; // indices of the upper-bounded variables

   // function evaluations
   [[nodiscard]] virtual double get_objective_multiplier() const = 0;
   [[nodiscard]] virtual double evaluate_objective(Iterate& iterate) const = 0;
   virtual void evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const = 0;
   virtual void evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const = 0;
   virtual void evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const = 0;
   virtual void evaluate_lagrangian_hessian(const std::vector<double>& x, const std::vector<double>& multipliers, SymmetricMatrix<double>& hessian) const = 0;

   [[nodiscard]] double compute_complementarity_error(const std::vector<double>& primals, const std::vector<double>& constraints,
         const std::vector<double>& constraint_multipliers, const std::vector<double>& lower_bounds_multipliers,
         const std::vector<double>& upper_bounds_multipliers) const;
   [[nodiscard]] double compute_dual_constraint_violation(const std::vector<double>& primals, const std::vector<double>& constraint_multipliers,
         const std::vector<double>& lower_bounds_multipliers, const std::vector<double>& upper_bounds_multipliers);

   [[nodiscard]] size_t get_number_original_variables() const;
   [[nodiscard]] virtual double get_variable_lower_bound(size_t i) const = 0;
   [[nodiscard]] virtual double get_variable_upper_bound(size_t i) const = 0;
   [[nodiscard]] virtual double get_constraint_lower_bound(size_t j) const = 0;
   [[nodiscard]] virtual double get_constraint_upper_bound(size_t j) const = 0;
   // relaxed bounds
   [[nodiscard]] double get_variable_lower_bound(size_t i, double relaxation_factor) const;
   [[nodiscard]] double get_variable_upper_bound(size_t i, double relaxation_factor) const;

   [[nodiscard]] virtual size_t get_maximum_number_objective_gradient_nonzeros() const = 0;
   [[nodiscard]] virtual size_t get_maximum_number_jacobian_nonzeros() const = 0;
   [[nodiscard]] virtual size_t get_maximum_number_hessian_nonzeros() const = 0;
};

inline NonlinearProblem::NonlinearProblem(const Model& model, size_t number_variables, size_t number_constraints):
      model(model), number_variables(number_variables), number_constraints(number_constraints) {
}

inline bool NonlinearProblem::is_constrained() const {
   return (0 < this->number_constraints);
}

inline size_t NonlinearProblem::get_number_original_variables() const {
   return this->model.number_variables;
}

// relaxed bounds
inline double NonlinearProblem::get_variable_lower_bound(size_t i, double relaxation_factor) const {
   const double lower_bound = this->get_variable_lower_bound(i);
   return lower_bound - relaxation_factor*std::max(1., std::abs(lower_bound));
}

inline double NonlinearProblem::get_variable_upper_bound(size_t i, double relaxation_factor) const {
   const double upper_bound = this->get_variable_upper_bound(i);
   return upper_bound + relaxation_factor * std::max(1., std::abs(upper_bound));
}

// complementary slackness error
inline double NonlinearProblem::compute_complementarity_error(const std::vector<double>& primals, const std::vector<double>& constraints,
      const std::vector<double>& constraint_multipliers, const std::vector<double>& lower_bounds_multipliers,
      const std::vector<double>& upper_bounds_multipliers) const {
   double error = 0.;
   // bound constraints
   for (size_t i: Range(this->number_variables)) {
      if (0. < lower_bounds_multipliers[i]) {
         //error += std::abs(lower_bounds_multipliers[i] * (primals[i] - this->get_variable_lower_bound(i)));
         //std::cout << lower_bounds_multipliers[i] << "*(" << primals[i] << " - " << this->get_variable_lower_bound(i) << ")\n";
         error = std::max(error, std::abs(lower_bounds_multipliers[i] * (primals[i] - this->get_variable_lower_bound(i))));
         //std::cout << "Error is " << error << '\n';
      }
      if (0. < upper_bounds_multipliers[i]) {
         //error += std::abs(upper_bounds_multipliers[i] * (primals[i] - this->get_variable_upper_bound(i)));
         //std::cout << upper_bounds_multipliers[i] << "*(" << primals[i] << " - " << this->get_variable_upper_bound(i) << ")\n";
         error = std::max(error, std::abs(upper_bounds_multipliers[i] * (primals[i] - this->get_variable_upper_bound(i))));
         //std::cout << "Error is " << error << '\n';
      }
   }
   //std::cout << "Error is " << error << '\n';
   // constraints
   this->inequality_constraints.for_each_index([&](size_t j) {
      if (0. < constraint_multipliers[j]) { // lower bound
         //error += std::abs(constraint_multipliers[j] * (constraints[j] - this->get_constraint_lower_bound(j)));
         error = std::max(error, std::abs(constraint_multipliers[j] * (constraints[j] - this->get_constraint_lower_bound(j))));
      }
      else if (constraint_multipliers[j] < 0.) { // upper bound
         //error += std::abs(constraint_multipliers[j] * (constraints[j] - this->get_constraint_upper_bound(j)));
         error = std::max(error, std::abs(constraint_multipliers[j] * (constraints[j] - this->get_constraint_upper_bound(j))));
      }
   });
   return error;
}

inline double NonlinearProblem::compute_dual_constraint_violation(const std::vector<double>& /*primals*/,
      const std::vector<double>& /*constraint_multipliers*/, const std::vector<double>& /*lower_bounds_multipliers*/,
            const std::vector<double>& /*upper_bounds_multipliers*/) {
   return 0.;
}

#endif // UNO_NONLINEARPROBLEM_H
