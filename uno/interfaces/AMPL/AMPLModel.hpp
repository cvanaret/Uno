#ifndef UNO_AMPLMODEL_H
#define UNO_AMPLMODEL_H

#include <vector>
#include "optimization/Model.hpp"

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

   // variables
   [[nodiscard]] double get_variable_lower_bound(size_t i) const override;
   [[nodiscard]] double get_variable_upper_bound(size_t i) const override;

   // objective
   [[nodiscard]] double evaluate_objective(const std::vector<double>& x) const override;
   void evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const override;

   // constraints
   void evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const override;
   void evaluate_constraint_gradient(const std::vector<double>& x, size_t j, SparseVector<double>& gradient) const override;
   void evaluate_constraint_jacobian(const std::vector<double>& x, std::vector<SparseVector<double>>& constraint_jacobian) const override;
   [[nodiscard]] double get_constraint_lower_bound(size_t j) const override;
   [[nodiscard]] double get_constraint_upper_bound(size_t j) const override;

   // Hessian
   void evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
         SymmetricMatrix& hessian) const override;

   [[nodiscard]] ConstraintType get_variable_status(size_t i) const override;
   [[nodiscard]] FunctionType get_constraint_type(size_t j) const override;
   [[nodiscard]] ConstraintType get_constraint_status(size_t j) const override;
   [[nodiscard]] size_t get_maximum_number_hessian_nonzeros() const override;

   void get_initial_primal_point(std::vector<double>& x) const override;
   void get_initial_dual_point(std::vector<double>& multipliers) const override;

private:
   // private constructor to pass the dimensions to the Problem base constructor
   AMPLModel(const std::string& file_name, ASL* asl);

   ASL* asl_; /*!< Instance of the AMPL Solver Library class */
   std::vector<Range> variables_bounds;
   std::vector<Range> constraint_bounds;
   std::vector<ConstraintType> variable_status; /*!< Status of the variables (EQUALITY, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES) */
   std::vector<FunctionType> constraint_type; /*!< Types of the constraints (LINEAR, QUADRATIC, NONLINEAR) */
   std::vector<ConstraintType> constraint_status; /*!< Status of the constraints (EQUAL_BOUNDS, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES,
 * UNBOUNDED) */

   std::vector<double> ampl_tmp_gradient{};
   std::vector<double> ampl_tmp_hessian{};

   void generate_variables();
   void generate_constraints();
   void set_function_types(std::string file_name);

   void initialize_lagrangian_hessian();
   [[nodiscard]] size_t compute_hessian_number_nonzeros(double objective_multiplier, const std::vector<double>& multipliers) const;
};

#endif // UNO_AMPLMODEL_H
