#ifndef AMPLMODEL_H
#define AMPLMODEL_H

#include <vector>
#include <map>
#include "Problem.hpp"
#include "Constraint.hpp"

extern "C" {
	#include "asl_pfgh.h"
	#include "getstub.h"
}

#define UNCERTAIN_SUFFIX "uncertain"
#define UNCERTAINTY_SET_SUFFIX "uncertainty_set"

/*! \class AMPLModel
* \brief AMPL model
*
*  Description of an AMPL model
*/
class AMPLModel: public Problem {
	public:
		/*!
         *  Constructor that takes the list of arguments of the command line
         */
		AMPLModel(std::string file_name);
		/*!
         *  Destructor
         */
		~AMPLModel();
		
		/* objective */
		double objective(std::vector<double> x);
		std::vector<double> objective_dense_gradient(std::vector<double> x);
		std::map<int,double> objective_sparse_gradient(std::vector<double> x);
		
		/* variables */
		std::vector<bool> variable_uncertain;
		
		/* constraints */
		std::vector<bool> constraint_is_uncertainty_set;
		double evaluate_constraint(int j, std::vector<double> x);
		std::vector<double> evaluate_constraints(std::vector<double> x);
		std::vector<double> constraint_dense_gradient(int j, std::vector<double> x);
		std::map<int,double> constraint_sparse_gradient(int j, std::vector<double> x);
		
		std::vector<int> jacobian_sparsity;
		std::vector<std::vector<double> > constraints_jacobian_dense(std::vector<double> x);
		void create_jacobian_sparsity();
		
		/* Hessian */
		Matrix lagrangian_hessian(std::vector<double> x, double objective_multiplier, std::vector<double> multipliers);
		
		std::vector<double> primal_initial_solution();
		std::vector<double> dual_initial_solution();
		
	private:
		ASL_pfgh* asl_; /*!< Instance of the AMPL Solver Library class */
		
		/*!
         *  Generate the variables
         */
		void generate_variables();
		
		/*!
         *  Generate the objective function
         */
		void initialize_objective();
		
		/*!
         *  Generate the constraints
         */
		void generate_constraints();
		
		/*!
         *  Determine the type (linear, quadratic, nonlinear) of the AMPL functions
         * 
         * \param file_name: list of arguments of the command line
         * \param option_info: AMPL options
         */
		void set_function_types(std::string file_name, Option_Info* option_info);
		
		/*!
         *  Generate the Hessian of the Lagrangian
         */
		void initialize_lagrangian_hessian();
};

#endif // AMPLMODEL_H
