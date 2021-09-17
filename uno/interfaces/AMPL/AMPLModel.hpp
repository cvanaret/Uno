#ifndef AMPLMODEL_H
#define AMPLMODEL_H

#include <vector>
#include "optimization/Problem.hpp"
#include "optimization/Constraint.hpp"

extern "C" {
#include "asl_pfgh.h"
#include "getstub.h"
}

/*! \class AMPLModel
 * \brief AMPL model
 *
 *  Description of an AMPL model
 */
class AMPLModel : public Problem {
public:
   explicit AMPLModel(std::string file_name);
   ~AMPLModel() override;

   /* objective */
   [[nodiscard]] double evaluate_objective(const std::vector<double>& x) const override;
   void evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const override;

   /* constraints */
   //std::vector<bool> constraint_is_uncertainty_set;
   [[nodiscard]] double evaluate_constraint(int j, const std::vector<double>& x) const override;
   void evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const override;
   void evaluate_constraint_gradient(const std::vector<double>& x, size_t j, SparseVector<double>& gradient) const override;
   void evaluate_constraint_jacobian(const std::vector<double>& x, std::vector<SparseVector<double>>& constraint_jacobian) const override;

   /* Hessian */
   void evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
         CSCSymmetricMatrix& hessian) const override;
   void evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
         COOSymmetricMatrix& hessian) const override;

   void set_initial_primal_point(std::vector<double>& x) override;
   void set_initial_dual_point(std::vector<double>& multipliers) override;

private:
   // private constructor to pass the dimensions to the Problem base constructor
   AMPLModel(std::string file_name, ASL* asl);

   ASL* asl_; /*!< Instance of the AMPL Solver Library class */
   std::vector<double> ampl_tmp_gradient{};
   std::vector<double> ampl_tmp_hessian{};

   void generate_variables();
   void generate_constraints();
   void set_function_types(std::string file_name);
   void initialize_lagrangian_hessian();
   size_t compute_hessian_number_nonzeros(double objective_multiplier, const std::vector<double>& multipliers) const;
   void generate_sparsity_pattern(CSCSymmetricMatrix& hessian, size_t number_non_zeros) const;
};

#endif // AMPLMODEL_H
