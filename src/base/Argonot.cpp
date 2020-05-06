#include <iostream>
#include <ctime>
#include <cmath>
#include "Argonot.hpp"
#include "Iterate.hpp"
#include "Logger.hpp"

Argonot::Argonot(GlobalizationMechanism& globalization_mechanism, int max_iterations): globalization_mechanism(globalization_mechanism), max_iterations(max_iterations) {
}

Result Argonot::solve(Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
    std::clock_t c_start = std::clock();
    int major_iterations = 0, minor_iterations = 0;

    /* use the current point to initialize the strategies and generate the initial point */
    Iterate current_iterate = this->globalization_mechanism.initialize(problem, x, multipliers);
    
    /* preprocessing phase: satisfy linear constraints */
    this->preprocessing(problem, current_iterate);
    
    INFO << "Problem " << problem.name << "\n";
    INFO << problem.number_variables << " variables, " << problem.number_constraints << " constraints\n";
    INFO << "Initial iterate\n" << current_iterate << "\n";

    try {
        /* check for convergence */
        while (!this->termination_criterion(current_iterate.status, major_iterations)) {
            major_iterations++;
            DEBUG << "\n\t\tARGONOT iteration " << major_iterations << "\n";
            INFO << "major: " << major_iterations << "\t";
            
            DEBUG << "Current point: "; print_vector(DEBUG, current_iterate.x);
            /* update the current point */
            current_iterate = this->globalization_mechanism.compute_acceptable_iterate(problem, current_iterate);
            minor_iterations += this->globalization_mechanism.number_iterations;
            INFO << "||c|| = " << current_iterate.constraint_residual << "\tf = " << current_iterate.objective << "\t";
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

    Result result = {problem.number_variables,
        problem.number_constraints,
        current_iterate,
        major_iterations,
        cpu_time,
        problem.number_eval_objective,
        problem.number_eval_constraints,
        problem.number_eval_jacobian,
        problem.number_eval_hessian,
        this->globalization_mechanism.globalization_strategy.subproblem.number_subproblems_solved};
    return result;
}

bool Argonot::termination_criterion(OptimalityStatus current_status, int iteration) {
    return current_status != NOT_OPTIMAL || this->max_iterations <= iteration;
}

double Argonot::compute_KKT_error(Problem& problem, Iterate& iterate, double objective_mutiplier, std::string norm_value) {
    std::vector<double> lagrangian_gradient = iterate.lagrangian_gradient(problem, objective_mutiplier, iterate.multipliers);
    return norm(lagrangian_gradient, norm_value);
}

void Argonot::preprocessing(Problem& problem, Iterate& iterate) {
    std::cout << "Preprocessing phase: nothing implemented\n";
}

void Result::display() {
    std::cout << "\n";
    std::cout << "ARGONOT v1: optimization summary\n";
    std::cout << "==============================\n";

    std::cout << "Status:\t\t\t\t";
    if (this->solution.status == KKT_POINT) {
        std::cout << "Converged with KKT point\n";
    }
    else if (this->solution.status == FJ_POINT) {
        std::cout << "Converged with FJ point\n";
    }
    else if (this->solution.status == FEASIBLE_SMALL_STEP) {
        std::cout << "Converged with feasible small step\n";
    }
    else if (this->solution.status == INFEASIBLE_SMALL_STEP) {
        std::cout << "Converged with infeasible small step\n";
    }
    else { // NOT_OPTIMAL
        std::cout << "Irregular termination\n";
    }

    std::cout << "Objective value:\t\t" << this->solution.objective << "\n";
    std::cout << "Constraint residual:\t\t" << this->solution.constraint_residual << "\n";
    std::cout << "KKT residual:\t\t\t" << this->solution.KKT_residual << "\n";
    std::cout << "Complementarity residual:\t" << this->solution.complementarity_residual << "\n";

    std::cout << "Feasibility measure:\t\t" << this->solution.feasibility_measure << "\n";
    std::cout << "Optimality measure:\t\t" << this->solution.optimality_measure << "\n";
    
    std::cout << "Primal solution:\t\t"; print_vector(std::cout, this->solution.x);
    std::cout << "Lower bound multipliers:\t"; print_vector(std::cout, this->solution.multipliers.lower_bounds);
    std::cout << "Upper bound multipliers:\t"; print_vector(std::cout, this->solution.multipliers.upper_bounds);
    std::cout << "Constraint multipliers:\t\t"; print_vector(std::cout, this->solution.multipliers.constraints);

    std::cout << "CPU time:\t\t\t" << this->cpu_time << "s\n";
    std::cout << "Iterations:\t\t\t" << this->iteration << "\n";
    std::cout << "Objective evaluations:\t\t" << this->objective_evaluations << "\n";
    std::cout << "Constraints evaluations:\t" << this->constraint_evaluations << "\n";
    std::cout << "Jacobian evaluations:\t\t" << this->jacobian_evaluations << "\n";
    std::cout << "Hessian evaluations:\t\t" << this->hessian_evaluations << "\n";
    std::cout << "Number of subproblems solved:\t" << this->number_subproblems_solved << "\n";
    return;
}
