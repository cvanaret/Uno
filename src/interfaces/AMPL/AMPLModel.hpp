#ifndef AMPLMODEL_H
#define AMPLMODEL_H

#include <vector>
#include <map>
#include <unordered_map>
#include "Problem.hpp"
#include "Constraint.hpp"

extern "C" {
#include "asl_pfgh.h"
#include "getstub.h"
}

//#define UNCERTAIN_SUFFIX "uncertain"
//#define UNCERTAINTY_SET_SUFFIX "uncertainty_set"

/*! \class AMPLModel
 * \brief AMPL model
 *
 *  Description of an AMPL model
 */
class AMPLModel : public Problem {
public:
    AMPLModel(std::string file_name, int fortran_indexing);
    ~AMPLModel();

    /* objective */
    double objective(const std::vector<double>& x) const override;
    std::vector<double> objective_dense_gradient(std::vector<double>& x) const override;
    SparseGradient objective_sparse_gradient(std::vector<double>& x) const override;
    
    /* constraints */
    //std::vector<bool> constraint_is_uncertainty_set;
    double evaluate_constraint(int j, std::vector<double>& x) const;
    std::vector<double> evaluate_constraints(std::vector<double>& x) const;
    std::vector<double> constraint_dense_gradient(int j, std::vector<double>& x) const;
    void constraint_sparse_gradient(std::vector<double>& x, int j, SparseGradient& gradient) const;
    std::vector<SparseGradient> constraints_sparse_jacobian(std::vector<double>& x) const;

    /* Hessian */
    CSCMatrix lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers) const override;
    //CSCMatrix lagrangian_hessian(std::vector<double>& x, double objective_multiplier, std::vector<double>& multipliers, std::vector<double>& hessian);

    std::vector<double> primal_initial_solution();
    std::vector<double> dual_initial_solution();

private:
    // private constructor to pass the dimensions to the Problem base constructor
    AMPLModel(std::string file_name, ASL_pfgh* asl, int fortran_indexing);

    ASL_pfgh* asl_; /*!< Instance of the AMPL Solver Library class */
    int fortran_indexing;
    std::vector<double> ampl_tmp_gradient_;

    void generate_variables_();
    void initialize_objective_();
    void generate_constraints_();
    //void create_objective_variables_(ograd* ampl_variables);
    //void create_constraint_variables_(int j, cgrad* ampl_variables);
    void set_function_types_(std::string file_name);
    void initialize_lagrangian_hessian_();
};

#endif // AMPLMODEL_H
