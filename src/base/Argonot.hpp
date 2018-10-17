#ifndef ARGONOT_H
#define ARGONOT_H

#include "Problem.hpp"
#include "GlobalizationMechanism.hpp"

struct Result {
    int number_variables;
    int number_constraints;
    Iterate solution;
    int iteration;
    double cpu_time;
    int objective_evaluations;
    int constraint_evaluations;
    int jacobian_evaluations;
    int hessian_evaluations;
    int number_subproblems_solved;

    void display();
};

/*! \class Argonot
* \brief Argonot
*
*  Argonot solver
*/
class Argonot {
	public:
		/*!
         *  Constructor
         * 
         * \param globalization_strategy: strategy to promote global convergence
         * \param tolerance: tolerance for termination criteria
         */
		Argonot(GlobalizationMechanism& globalization_mechanism, int max_iterations = 1000);
		
		/*!
         *  Solve a given problem with initial primal and dual variables
         * 
         * \param problem: optimization problem
         * \param x: primal variables
         * \param multipliers: Lagrange multipliers/dual variables
         */
		Result solve(Problem& problem, std::vector<double>& x, std::vector<double>& bound_multipliers, std::vector<double>& constraint_multipliers);
		
		GlobalizationMechanism& globalization_mechanism; /*!< Step control strategy (trust region or line-search) */
		int max_iterations; /*!< Maximum number of iterations */
	
	private:
		/*!
         *  Determine whether the optimization process is over
         * 
         * \param is_optimal: optimality status
         * \param iteration: current iteration number
         */
		bool termination_criterion(OptimalityStatus is_optimal, int iteration);
		
		/*!
         *  Determine if the current point satisfies optimality criteria
         * 
         * \param problem: optimization problem
         * \param phase: current phase (optimality or feasibility restoration)
         * \param current_iterate: current point and its evaluations
         */
		OptimalityStatus optimality_test(Problem& problem, Phase& phase, Iterate& current_iterate);
		
		/*!
         *  Compute the KKT error
         * 
         * \param problem: optimization problem
         * \param phase: current phase (optimality or feasibility restoration)
         * \param current_iterate: current point and its evaluations

         */
		double compute_KKT_error(Problem& problem, Phase& phase, Iterate& current_iterate);
			
		/*!
         *  Compute the complementarity error
         * 
         * \param problem: optimization problem
         * \param phase: current phase (optimality or feasibility restoration)
         * \param current_iterate: current point and its evaluations
         */
		double compute_complementarity_error(const Problem& problem, Phase& phase, Iterate& current_iterate);
};

#endif // ARGONOT_H
