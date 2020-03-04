#include <iostream>
#include <ctime>
#include <cmath>
#include "Argonot.hpp"
#include "Iterate.hpp"
#include "Logger.hpp"

Argonot::Argonot(GlobalizationMechanism& globalization_mechanism, int max_iterations) :
globalization_mechanism(globalization_mechanism), max_iterations(max_iterations) {
}

Result Argonot::solve(Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
    std::clock_t c_start = std::clock();

    int major_iterations = 0, minor_iterations = 0;

    INFO << "Problem " << problem.name << "\n";
    INFO << problem.number_variables << " variables, " << problem.number_constraints << " constraints\n";

    /* use the current point to initialize the strategies and generate the initial point */
    Iterate current_iterate = this->globalization_mechanism.initialize(problem, x, multipliers);
    INFO << "Initial iterate\n" << current_iterate << "\n";

    try {
        /* check for convergence */
        while (!this->termination_criterion(current_iterate.status, major_iterations)) {
            major_iterations++;
            DEBUG << "\n\t\tARGONOT iteration " << major_iterations << "\n";
            INFO << "major: " << major_iterations << "\t";

            /* update the current point */

            current_iterate = this->globalization_mechanism.compute_iterate(problem, current_iterate);
            minor_iterations += this->globalization_mechanism.number_iterations;

            INFO << "constraints: " << current_iterate.residual << "\tobjective: " << current_iterate.objective << "\t";
            INFO << "status: " << current_iterate.status << "\n";
            DEBUG << "Next iterate\n" << current_iterate;
        }
    }    catch (std::out_of_range& exception) {
        ERROR << exception.what();
    }    catch (std::invalid_argument& exception) {
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

std::vector<double> Argonot::compute_lagrangian_gradient(Problem& problem, Iterate& current_iterate, double objective_multiplier, Multipliers& multipliers) {
    std::vector<double> lagrangian_gradient(problem.number_variables);

    /* objective gradient */
    if (objective_multiplier != 0.) {
        current_iterate.compute_objective_gradient(problem);
        
        /* scale the objective gradient */
        for (std::pair<int, double> term : current_iterate.objective_gradient) {
            int variable_index = term.first;
            double derivative = term.second;
            lagrangian_gradient[variable_index] += objective_multiplier*derivative;
        }
    }
    /* bound constraints */
    for (int i = 0; i < problem.number_variables; i++) {
        lagrangian_gradient[i] += - multipliers.lower_bounds[i] - multipliers.upper_bounds[i];
    }
    /* constraints */
    current_iterate.compute_constraints_jacobian(problem);

    for (int j = 0; j < problem.number_constraints; j++) {
        double multiplier_j = multipliers.constraints[j];
        if (multiplier_j != 0.) {
            for (std::pair<int, double> term : current_iterate.constraints_jacobian[j]) {
                int variable_index = term.first;
                double derivative = term.second;
                lagrangian_gradient[variable_index] -= multiplier_j*derivative;
            }
        }
    }
    return lagrangian_gradient;
}

double Argonot::compute_KKT_error(Problem& problem, Iterate& iterate, double objective_mutiplier, std::string norm) {
    std::vector<double> lagrangian_gradient = Argonot::compute_lagrangian_gradient(problem, iterate, objective_mutiplier, iterate.multipliers);
    if (norm == "inf") {
        return norm_inf(lagrangian_gradient);
    }
    else if (norm == "l2") {
        return norm_2(lagrangian_gradient);
    }
    else if (norm == "l1") {
        return norm_1(lagrangian_gradient);
    }
    else {
        throw std::out_of_range("The norm is not known");
    }
}

/* complementary slackness error. Use abs/1e-8 to safeguard */
double Argonot::compute_complementarity_error(const Problem& problem, Iterate& iterate) {
    double complementarity_error = 0.;
	
    /* bound constraints */
    for (int i = 0; i < problem.number_variables; i++) {
        if (-INFINITY < problem.variables_bounds[i].lb) {
            complementarity_error += std::abs(iterate.multipliers.lower_bounds[i] * (iterate.x[i] - problem.variables_bounds[i].lb));
        }
        else if (problem.variables_bounds[i].ub < INFINITY) {
            complementarity_error += std::abs(iterate.multipliers.upper_bounds[i] * (iterate.x[i] - problem.variables_bounds[i].ub));
        }
    }
    /* constraints */
    for (int j = 0; j < problem.number_constraints; j++) {
        double multiplier_j = iterate.multipliers.constraints[j];

        if (-INFINITY < problem.constraints_bounds[j].lb) {
            complementarity_error += std::abs(multiplier_j * (iterate.constraints[j] - problem.constraints_bounds[j].lb));
        }
        else if (problem.constraints_bounds[j].ub < INFINITY) {
            complementarity_error += std::abs(multiplier_j * (iterate.constraints[j] - problem.constraints_bounds[j].ub));
        }
    }
    return complementarity_error;
}

void Result::display() {
    std::cout << "\n";
    std::cout << "ARGONOT v1: optimization summary\n";
    std::cout << "==============================\n";

    std::cout << "Status:\t\t\t";
    if (this->solution.status == KKT_POINT) {
        std::cout << "Feasible KKT point found\n";
    } else if (this->solution.status == FJ_POINT) {
        std::cout << "Infeasible FJ point found\n";
    } else if (this->solution.status == FEASIBLE_SMALL_STEP) {
        std::cout << "Feasible small step\n";
    } else if (this->solution.status == INFEASIBLE_SMALL_STEP) {
        std::cout << "Infeasible small step\n";
    } else { // NOT_OPTIMAL
        std::cout << "Irregular termination\n";
    }

    std::cout << "Objective value:\t" << this->solution.objective << "\n";
    std::cout << "Constraint residual:\t" << this->solution.residual << "\n";
    std::cout << "KKT error:\t\t" << this->solution.KKTerror << "\n";
    std::cout << "Complementarity error:\t" << this->solution.complementarity_error << "\n";

    std::cout << "Primal solution:\t";
    print_vector(std::cout, this->solution.x);

    std::cout << "Lower bound multipliers:\t";
    print_vector(std::cout, this->solution.multipliers.lower_bounds);
    std::cout << "Upper bound multipliers:\t";
    print_vector(std::cout, this->solution.multipliers.upper_bounds);
    std::cout << "Constraint multipliers:\t";
    print_vector(std::cout, this->solution.multipliers.constraints);

    std::cout << "CPU time:\t\t" << this->cpu_time << "s\n";
    std::cout << "Iterations:\t\t" << this->iteration << "\n";
    std::cout << "Objective evaluations:\t\t" << this->objective_evaluations << "\n";
    std::cout << "Constraints evaluations:\t\t" << this->constraint_evaluations << "\n";
    std::cout << "Jacobian evaluations:\t\t" << this->jacobian_evaluations << "\n";
    std::cout << "Hessian evaluations:\t\t" << this->hessian_evaluations << "\n";
    std::cout << "Number of subproblems solved:\t\t" << this->number_subproblems_solved << "\n";
    return;
}
