// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_AMPLMODEL_H
#define UNO_AMPLMODEL_H

#include <vector>
#include "optimization/Model.hpp"
#include "linear_algebra/RectangularMatrix.hpp"

// include AMPL Solver Library (ASL)
extern "C" {
#include "asl_pfgh.h"
#include "getstub.h"
}

/*! \class AMPLModel
 * \brief AMPL model
 *
 *  Description of an AMPL model
 */
class AMPLModel: public Model {
public:
   explicit AMPLModel(const std::string& file_name);
   ~AMPLModel() override;

   // objective
   [[nodiscard]] double evaluate_objective(const std::vector<double>& x) const override;
   void evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const override;
   // constraints
   void evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const override;
   void evaluate_constraint_gradient(const std::vector<double>& x, size_t j, SparseVector<double>& gradient) const override;
   void evaluate_constraint_jacobian(const std::vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const override;
   // Hessian
   void evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
         SymmetricMatrix<double>& hessian) const override;

   [[nodiscard]] double get_variable_lower_bound(size_t i) const override;
   [[nodiscard]] double get_variable_upper_bound(size_t i) const override;
   [[nodiscard]] BoundType get_variable_bound_type(size_t i) const override;
   [[nodiscard]] double get_constraint_lower_bound(size_t j) const override;
   [[nodiscard]] double get_constraint_upper_bound(size_t j) const override;
   [[nodiscard]] FunctionType get_constraint_type(size_t j) const override;
   [[nodiscard]] BoundType get_constraint_bound_type(size_t j) const override;

   [[nodiscard]] size_t get_number_objective_gradient_nonzeros() const override;
   [[nodiscard]] size_t get_number_jacobian_nonzeros() const override;
   [[nodiscard]] size_t get_number_hessian_nonzeros() const override;

   void get_initial_primal_point(std::vector<double>& x) const override;
   void get_initial_dual_point(std::vector<double>& multipliers) const override;
   void postprocess_solution(Iterate& iterate, TerminationStatus termination_status) const override;

   [[nodiscard]] const std::vector<size_t>& get_linear_constraints() const override;

private:
   // private constructor to pass the dimensions to the Model base constructor
   AMPLModel(const std::string& file_name, ASL* asl);

   // mutable: can be modified by const methods (internal state not seen by user)
   mutable ASL* asl; /*!< Instance of the AMPL Solver Library class */
   mutable std::vector<double> ampl_tmp_gradient{};
   mutable std::vector<double> ampl_tmp_hessian{};

   std::vector<Interval> variables_bounds;
   std::vector<Interval> constraint_bounds;
   std::vector<BoundType> variable_status; /*!< Status of the variables (EQUALITY, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES) */
   std::vector<FunctionType> constraint_type; /*!< Types of the constraints (LINEAR, QUADRATIC, NONLINEAR) */
   std::vector<BoundType> constraint_status; /*!< Status of the constraints (EQUAL_BOUNDS, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES,
 * UNBOUNDED) */

   std::vector<size_t> linear_constraints;

   void generate_variables();
   void generate_constraints();

   void set_number_hessian_nonzeros();
   [[nodiscard]] size_t compute_hessian_number_nonzeros(double objective_multiplier, const std::vector<double>& multipliers) const;
};

#endif // UNO_AMPLMODEL_H
