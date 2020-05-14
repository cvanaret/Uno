#include <iostream>
#include <ctime>
#include <cmath>
#include "Argonot.hpp"
#include "Iterate.hpp"
#include "Logger.hpp"
#include "BQPDSolver.hpp"

Argonot::Argonot(GlobalizationMechanism& globalization_mechanism, int max_iterations): globalization_mechanism(globalization_mechanism), max_iterations(max_iterations) {
}

Result Argonot::solve(Problem& problem, std::vector<double>& x, Multipliers& multipliers, bool preprocessing) {
    std::clock_t c_start = std::clock();
    int major_iterations = 0;
    
    INFO << "Problem " << problem.name << "\n";
    INFO << problem.number_variables << " variables, " << problem.number_constraints << " constraints\n";
    
    /* project x into the bounds */
    Subproblem::project_point_in_bounds(x, problem.variables_bounds);
    
    if (preprocessing) {
        /* preprocessing phase: satisfy linear constraints */
        this->preprocessing(problem, x, multipliers);
    }
    
    /* use the current point to initialize the strategies and generate the initial point */
    Iterate current_iterate = this->globalization_mechanism.initialize(problem, x, multipliers);
    DEBUG << "Initial iterate\n" << current_iterate << "\n";
    
    try {
        /* check for convergence */
        while (!this->termination_criterion(current_iterate.status, major_iterations)) {
            major_iterations++;
            DEBUG << "\n\t\tARGONOT iteration " << major_iterations << "\n";
            INFO << "major: " << major_iterations << "\t";
            
            DEBUG << "Current point: "; print_vector(DEBUG, current_iterate.x);
            /* update the current point */
            current_iterate = this->globalization_mechanism.compute_acceptable_iterate(problem, current_iterate);
            INFO << "||c|| = " << current_iterate.residuals.constraints << "\tf = " << current_iterate.objective << "\t";
            INFO << "η = " << current_iterate.feasibility_measure << "\tω = " << current_iterate.optimality_measure << "\n";

            DEBUG << "Next iterate\n" << current_iterate;
        }
    }
    catch (std::invalid_argument& exception) {
        ERROR << exception.what();
    }
    catch (std::runtime_error& exception) {
        ERROR << exception.what();
    }
    std::clock_t c_end = std::clock();
    double cpu_time = (c_end - c_start) / (double) CLOCKS_PER_SEC;
    
    int number_subproblems_solved = this->globalization_mechanism.globalization_strategy.subproblem.number_subproblems_solved;
    Result result = {current_iterate,
        problem.number_variables,
        problem.number_constraints,
        major_iterations,
        cpu_time,
        problem.number_eval_objective,
        problem.number_eval_constraints,
        problem.number_eval_jacobian,
        problem.number_eval_hessian,
        number_subproblems_solved};
    return result;
}

bool Argonot::termination_criterion(OptimalityStatus current_status, int iteration) {
    return current_status != NOT_OPTIMAL || this->max_iterations <= iteration;
}

void Argonot::preprocessing(Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
    INFO << "Preprocessing phase:\n";
    
    /* linear constraints */
    INFO << "The problem has " << problem.linear_constraints.size() << " linear constraints\n";
    if (0 < problem.linear_constraints.size()) {
        std::vector<double> constraints = problem.evaluate_constraints(x);
        
        int infeasible_linear_constraints = 0;
        for (std::pair<int, int> element: problem.linear_constraints) {
            int j = element.first;
            if (constraints[j] < problem.constraint_bounds[j].lb || problem.constraint_bounds[j].ub < constraints[j]) {
                infeasible_linear_constraints++;
            } 
        }
        INFO << "There are " << infeasible_linear_constraints << " infeasible linear constraints at the initial point\n";
        
        if (0 < infeasible_linear_constraints) {
            INFO << "Current point: "; print_vector(INFO, x);
            int number_constraints = (int) problem.linear_constraints.size();
            BQPDSolver solver(problem.number_variables, number_constraints, problem.number_variables);
            
            int fortran_indexing = 1;
            CSCMatrix hessian = CSCMatrix::identity(problem.number_variables, fortran_indexing);
            std::map<int, double> linear_objective; // empty
            std::vector<double> d0(problem.number_variables);
            // constraints Jacobian
            std::vector<std::map<int, double> > constraints_jacobian(number_constraints);
            for (std::pair<int, int> element: problem.linear_constraints) {
                int j = element.first;
                int linear_constraint_index = element.second;
                problem.constraint_sparse_gradient(x, j, constraints_jacobian[linear_constraint_index]);
            }
            // variables bounds
            std::vector<Range> variables_bounds(problem.number_variables);
            for (int i = 0; i < problem.number_variables; i++) {
                variables_bounds[i] = {problem.variables_bounds[i].lb - x[i], problem.variables_bounds[i].ub - x[i]};
            }
            // constraints bounds
            std::vector<Range> constraints_bounds(number_constraints);
            for (std::pair<int, int> element: problem.linear_constraints) {
                int j = element.first;
                int linear_constraint_index = element.second;
                constraints_bounds[linear_constraint_index] = {problem.constraint_bounds[j].lb - constraints[j], problem.constraint_bounds[j].ub - constraints[j]};
            }
            SubproblemSolution solution = solver.solve_QP(variables_bounds, constraints_bounds, linear_objective, constraints_jacobian, hessian, d0);
            if (solution.status == INFEASIBLE) {
                throw std::runtime_error("Linear constraints cannot be satisfied");
            }
            
            std::vector<double> feasible_x = add_vectors(x, solution.x, 1.);
            x = feasible_x;
            // copy bound multipliers
            multipliers.lower_bounds = solution.multipliers.lower_bounds;
            multipliers.upper_bounds = solution.multipliers.upper_bounds;
            // copy constraint multipliers
            for (std::pair<int, int> element: problem.linear_constraints) {
                int j = element.first;
                int linear_constraint_index = element.second;
                multipliers.constraints[j] = solution.multipliers.constraints[linear_constraint_index];
            }
            INFO << "Linear feasible initial point: "; print_vector(INFO, x);
            INFO << "\n";
        }
    }
    return;
}

void Result::display(bool print_solution) {
    INFO << "\n";
    INFO << "ARGONOT v1: optimization summary\n";
    INFO << "==============================\n";

    INFO << "Status:\t\t\t\t";
    if (this->solution.status == KKT_POINT) {
        INFO << "Converged with KKT point\n";
    }
    else if (this->solution.status == FJ_POINT) {
        INFO << "Converged with FJ point\n";
    }
    else if (this->solution.status == FEASIBLE_SMALL_STEP) {
        INFO << "Converged with feasible small step\n";
    }
    else if (this->solution.status == INFEASIBLE_SMALL_STEP) {
        INFO << "Converged with infeasible small step\n";
    }
    else { // NOT_OPTIMAL
        INFO << "Irregular termination\n";
    }

    INFO << "Objective value:\t\t" << this->solution.objective << "\n";
    INFO << "Constraint residual:\t\t" << this->solution.residuals.constraints << "\n";
    INFO << "KKT residual:\t\t\t" << this->solution.residuals.KKT << "\n";
    INFO << "FJ residual:\t\t\t" << this->solution.residuals.FJ << "\n";
    INFO << "Complementarity residual:\t" << this->solution.residuals.complementarity << "\n";

    INFO << "Feasibility measure:\t\t" << this->solution.feasibility_measure << "\n";
    INFO << "Optimality measure:\t\t" << this->solution.optimality_measure << "\n";
    
    if (print_solution) {
        INFO << "Primal solution:\t\t"; print_vector(INFO, this->solution.x);
        INFO << "Lower bound multipliers:\t"; print_vector(INFO, this->solution.multipliers.lower_bounds);
        INFO << "Upper bound multipliers:\t"; print_vector(INFO, this->solution.multipliers.upper_bounds);
        INFO << "Constraint multipliers:\t\t"; print_vector(INFO, this->solution.multipliers.constraints);
    }

    INFO << "CPU time:\t\t\t" << this->cpu_time << "s\n";
    INFO << "Iterations:\t\t\t" << this->iteration << "\n";
    INFO << "Objective evaluations:\t\t" << this->objective_evaluations << "\n";
    INFO << "Constraints evaluations:\t" << this->constraint_evaluations << "\n";
    INFO << "Jacobian evaluations:\t\t" << this->jacobian_evaluations << "\n";
    INFO << "Hessian evaluations:\t\t" << this->hessian_evaluations << "\n";
    INFO << "Number of subproblems solved:\t" << this->number_subproblems_solved << "\n";
    return;
}
