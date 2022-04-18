#ifndef UNO_NONLINEARPROBLEM_H
#define UNO_NONLINEARPROBLEM_H

#include <string>
#include <vector>
#include <map>
#include "optimization/Constraint.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/Model.hpp"
#include "linear_algebra/CSCSymmetricMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"

class NonlinearProblem {
public:
   NonlinearProblem(const Model& model, size_t number_variables, size_t number_constraints);
   virtual ~NonlinearProblem() = default;

   const Model& model;
   const size_t number_variables; /*!< Number of variables */
   const size_t number_constraints; /*!< Number of constraints */

   SparseVector<size_t> equality_constraints; /*!< inequality constraints */
   SparseVector<size_t> inequality_constraints; /*!< inequality constraints */
   SparseVector<size_t> linear_constraints;
   // lists of bounded variables
   std::vector<size_t> lower_bounded_variables{}; // indices of the lower-bounded variables
   std::vector<size_t> upper_bounded_variables{}; // indices of the upper-bounded variables

   // function evaluations
   [[nodiscard]] virtual double evaluate_objective(Iterate& iterate) const = 0;
   virtual void evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const = 0;
   virtual void evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const = 0;
   virtual void evaluate_constraint_jacobian(Iterate& iterate, std::vector<SparseVector<double>>& constraint_jacobian) const = 0;
   virtual void evaluate_lagrangian_hessian(const std::vector<double>& x, const std::vector<double>& multipliers, SymmetricMatrix& hessian) const = 0;
   void evaluate_lagrangian_gradient(Iterate& iterate, std::vector<double>& lagrangian_gradient) const;

   [[nodiscard]] size_t get_number_original_variables() const;
   [[nodiscard]] virtual double get_variable_lower_bound(size_t i) const = 0;
   [[nodiscard]] virtual double get_variable_upper_bound(size_t i) const = 0;
   [[nodiscard]] virtual double get_constraint_lower_bound(size_t j) const = 0;
   [[nodiscard]] virtual double get_constraint_upper_bound(size_t j) const = 0;

   [[nodiscard]] virtual size_t get_hessian_maximum_number_nonzeros() const = 0;

   [[nodiscard]] double compute_constraint_lower_bound_violation(double constraint, size_t j) const;
   [[nodiscard]] double compute_constraint_upper_bound_violation(double constraint, size_t j) const;
   [[nodiscard]] virtual double compute_constraint_violation(double constraint, size_t j) const;
   [[nodiscard]] double compute_constraint_violation(const std::vector<double>& constraints, Norm residual_norm) const;
   [[nodiscard]] double compute_constraint_violation(const std::vector<double>& constraints, const std::vector<size_t>& constraint_set,
         Norm residual_norm) const;

protected:
   size_t hessian_maximum_number_nonzeros{0}; /*!< Number of nonzero elements in the Hessian */
};

#endif // UNO_NONLINEARPROBLEM_H
