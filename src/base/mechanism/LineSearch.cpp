#include <cmath>
#include <assert.h>
#include "LineSearch.hpp"
#include "Logger.hpp"
#include "InteriorPoint.hpp"

LineSearch::LineSearch(GlobalizationStrategy& globalization_strategy, double tolerance, int max_iterations, double backtracking_ratio):
GlobalizationMechanism(globalization_strategy, tolerance, max_iterations), backtracking_ratio(backtracking_ratio), min_step_length(1e-9) {
}

Iterate LineSearch::initialize(Statistics& statistics, Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
    statistics.add_column("minor", Statistics::int_width, 2);
    statistics.add_column("LS step length", Statistics::double_width, 30);
    statistics.add_column("step norm", Statistics::double_width, 31);
    return this->globalization_strategy.initialize(statistics, problem, x, multipliers);
}

Iterate LineSearch::compute_acceptable_iterate(Statistics& statistics, Problem& problem, Iterate& current_iterate) {
    bool line_search_termination = false;
    /* compute the step */
    std::vector<Direction> directions = this->globalization_strategy.subproblem.compute_directions(problem, current_iterate);

    while (!line_search_termination) {
        /* step length follows the following sequence: 1, ratio, ratio^2, ratio^3, ... */
        this->step_length = 1.;
        bool is_accepted = false;
        this->number_iterations = 0;

        while (!this->termination_(is_accepted)) {
            assert (0 < this->step_length && this->step_length <= 1);
            this->number_iterations++;
            this->print_iteration_();

            /* check whether the trial step is accepted */
            std::optional<std::pair<Iterate, Direction> > acceptance_check = this->find_first_acceptable_direction_(statistics, problem, current_iterate, directions, this->step_length);
            if (acceptance_check.has_value()) {
                is_accepted = true;
                current_iterate = acceptance_check.value().first;
                Direction direction = acceptance_check.value().second;
                statistics.add_statistic("minor", this->number_iterations);
                statistics.add_statistic("LS step length", this->step_length);
                statistics.add_statistic("step norm", this->step_length*direction.norm);
            }
            if (!is_accepted) {
                /* decrease the step length */
                this->step_length *= this->backtracking_ratio;
            }
        }
        // if step length is too small, run restoration phase
        if (!is_accepted && this->step_length < this->min_step_length) {
            if (0. < current_iterate.feasibility_measure && directions[0].phase == OPTIMALITY) {
                // reset the line search with the restoration solution
                DEBUG << "Enter restoration feasibility phase\n";
                directions = this->globalization_strategy.subproblem.restore_feasibility(problem, current_iterate, directions[0]);
                this->step_length = 1.;
            }
            else {
                throw std::runtime_error("Line-search iteration limit reached");
            }
        }
        else {
            line_search_termination = true;
        }
    }
    return current_iterate;
}

bool LineSearch::termination_(bool is_accepted) {
    if (is_accepted) {
        return true;
    }
    else if (this->max_iterations < this->number_iterations) {
        return true;
    }
    return false;
}

void LineSearch::print_iteration_() {
    DEBUG << "\tLINE SEARCH iteration " << this->number_iterations << ", step_length " << this->step_length << "\n";
    return;
}

void LineSearch::print_acceptance_() {
    DEBUG << CYAN "LS trial point accepted\n" RESET;
    return;
}

void LineSearch::print_warning_(const char* message) {
    WARNING << RED << message << RESET "\n";
    return;
}

/*
 * Interpolation functions
 */

//double LineSearch::quadratic_interpolation(Problem& problem, Iterate& current_iterate, std::vector<double> direction, double step_length) {
//    std::cout << "Current point: ";
//    print_vector(std::cout, current_iterate.x);
//    std::cout << "Direction: ";
//    print_vector(std::cout, direction);
//
//    /* compute trial point */
//    std::vector<double> trial_point(problem.number_variables);
//    for (int i = 0; i < problem.number_variables; i++) {
//        trial_point[i] = current_iterate.x[i] + step_length * direction[i];
//    }
//    /* evaluate trial point */
//    double phi_alpha0 = problem.objective(trial_point);
//    std::cout << "phi(alpha0) = f(x + alpha0*p) = " << phi_alpha0 << "\n";
//    /* compute dot product */
//    current_iterate.compute_objective_gradient(problem);
//    double phi_prime_0 = dot(direction, current_iterate.objective_gradient);
//    std::cout << "phi'(0) = nabla f(x)^T p = " << phi_prime_0 << "\n";
//
//    /* compute the minimum of the quadratic */
//    double a = (phi_alpha0 - current_iterate.objective - phi_prime_0 * step_length) / (step_length * step_length);
//    double b = phi_prime_0;
//    std::cout << "a = " << a << ", b = " << b << "\n";
//    return this->minimize_quadratic(a, b);
//}
//
//double LineSearch::cubic_interpolation(Problem& problem, Iterate& current_iterate, std::vector<double> direction, double steplength1, double steplength2) {
//    std::cout << "Current point: ";
//    print_vector(std::cout, current_iterate.x);
//    std::cout << "Direction: ";
//    print_vector(std::cout, direction);
//
//    /* compute trial points */
//    std::vector<double> trial_point1(problem.number_variables);
//    std::vector<double> trial_point2(problem.number_variables);
//    for (int i = 0; i < problem.number_variables; i++) {
//        trial_point1[i] = current_iterate.x[i] + steplength1 * direction[i];
//        trial_point2[i] = current_iterate.x[i] + steplength2 * direction[i];
//    }
//    /* evaluate trial points */
//    double phi_alpha1 = problem.objective(trial_point1);
//    double phi_alpha2 = problem.objective(trial_point2);
//    std::cout << "phi(alpha1) = f(x + alpha1*p) = " << phi_alpha1 << "\n";
//    std::cout << "phi(alpha2) = f(x + alpha2*p) = " << phi_alpha2 << "\n";
//    /* compute dot product */
//    current_iterate.compute_objective_gradient(problem);
//    double phi_prime_0 = dot(direction, current_iterate.objective_gradient);
//    std::cout << "phi'(0) = nabla f(x)^T p = " << phi_prime_0 << "\n";
//
//    /* compute the minimum of the cubic */
//    double det = steplength1 * steplength1 * steplength2 * steplength2 * (steplength1 - steplength2);
//    std::cout << "Det = " << det << "\n";
//    double K1 = phi_alpha1 - current_iterate.objective - steplength1*phi_prime_0;
//    double K2 = phi_alpha2 - current_iterate.objective - steplength2*phi_prime_0;
//
//    double a = (steplength2 * steplength2 * K1 - steplength1 * steplength1 * K2) / det;
//    double b = (-steplength2 * steplength2 * steplength2 * K1 + steplength1 * steplength1 * steplength1 * K2) / det; //phi_prime_0;
//    double c = phi_prime_0;
//    std::cout << "a = " << a << ", b = " << b << ", c = " << c << "\n";
//    if (a == 0.) {
//        return this->minimize_quadratic(b, c);
//    }
//    else {
//        return this->minimize_cubic(a, b, c);
//    }
//}
//
///* return the minimum of x -> ax^2 + bx + R */
//double LineSearch::minimize_quadratic(double a, double b) {
//    return -b / (2. * a);
//}
//
///* return the minimum of x -> ax^3 + bx^2 + cx + R */
//double LineSearch::minimize_cubic(double a, double b, double c) {
//    return (-b + std::sqrt(b * b - 3. * a * c)) / (3. * a);
//}
