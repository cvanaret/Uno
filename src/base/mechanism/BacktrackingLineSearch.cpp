#include <cmath>
#include <cassert>
#include "BacktrackingLineSearch.hpp"
#include "Logger.hpp"
#include "InteriorPoint.hpp"

BacktrackingLineSearch::BacktrackingLineSearch(ConstraintRelaxationStrategy& constraint_relaxation_strategy, int max_iterations, double backtracking_ratio):
   GlobalizationMechanism(constraint_relaxation_strategy, max_iterations), step_length(1.), backtracking_ratio(backtracking_ratio),
   min_step_length(1e-9) {
}

Iterate BacktrackingLineSearch::initialize(Statistics& statistics, const Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
   statistics.add_column("LS step length", Statistics::double_width, 30);
   // generate the initial point
   Iterate first_iterate = this->relaxation_strategy.initialize(statistics, problem, x, multipliers);

   // preallocate trial_iterate
   this->trial_primals_.resize(first_iterate.x.size());

   return first_iterate;
}

std::pair<Iterate, Direction> BacktrackingLineSearch::compute_acceptable_iterate(Statistics& statistics, const Problem& problem, Iterate& current_iterate) {
   /* compute the directions */
   this->relaxation_strategy.subproblem.generate(problem, current_iterate, problem.objective_sign, INFINITY);
   Direction direction = this->relaxation_strategy.compute_feasible_direction(problem, current_iterate, INFINITY);

   bool line_search_termination = false;
   while (!line_search_termination) {
      /* step length follows the following sequence: 1, ratio, ratio^2, ratio^3, ... */
      this->step_length = 1.;
      this->number_iterations = 0;
      while (!this->termination_()) {
         assert (0 < this->step_length && this->step_length <= 1);
         this->number_iterations++;
         this->print_iteration_();

         // assemble the trial iterate
         Iterate trial_iterate = this->assemble_trial_iterate(problem, current_iterate, direction, this->step_length);

         // check whether the trial step is accepted
         if (this->relaxation_strategy.is_acceptable(statistics, problem, current_iterate, trial_iterate, direction, this->step_length)) {
            // compute the residuals
            trial_iterate.compute_objective(problem);
            this->relaxation_strategy.subproblem.compute_residuals(problem, trial_iterate, trial_iterate.multipliers, direction.objective_multiplier);

            this->add_statistics(statistics, direction);
            return std::make_pair(trial_iterate, direction);
         }
         /* decrease the step length */
         this->update_step_length();
      }
      // if step length is too small, run restoration phase
      if (this->step_length < this->min_step_length) {
         if (0. < current_iterate.progress.feasibility && !direction.is_relaxed) {
            // reset the line search with the restoration solution
            DEBUG << "Enter restoration feasibility phase\n";
            direction = this->relaxation_strategy.subproblem.restore_feasibility(problem, current_iterate, direction, INFINITY);
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
   throw std::runtime_error("Trust-region failed with an unexpected error");
}

void BacktrackingLineSearch::add_statistics(Statistics& statistics, const Direction& direction) {
   statistics.add_statistic("minor", this->number_iterations);
   statistics.add_statistic("LS step length", this->step_length);
   statistics.add_statistic("step norm", this->step_length * direction.norm);
}

void BacktrackingLineSearch::update_step_length() {
   this->step_length *= this->backtracking_ratio;
}

bool BacktrackingLineSearch::termination_() {
   return (this->max_iterations < this->number_iterations);
}

void BacktrackingLineSearch::print_iteration_() {
   DEBUG << "\tLINE SEARCH iteration " << this->number_iterations << ", step_length " << this->step_length << "\n";
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
