#ifndef IPM_H
#define IPM_H

#include <exception>
#include "Subproblem.hpp"
#include "solvers/linear/LinearSolver.hpp"
#include "HessianModel.hpp"

struct InteriorPointParameters {
   double tau_min;
   double k_sigma;
   double smax;
   double k_mu;
   double theta_mu;
   double k_epsilon;
};

struct UnstableInertiaCorrection : public std::exception {

   [[nodiscard]] const char* what() const throw() override {
      return "The inertia correction got unstable (delta_w > 1e40)";
   }
};

template<typename LinearSolverType>
class InteriorPoint : public Subproblem {
public:
   InteriorPoint(const Problem& problem, size_t max_number_variables, size_t number_constraints, const std::string& hessian_evaluation_method,
         double initial_barrier_parameter, double default_multiplier, double tolerance, bool use_trust_region);
   ~InteriorPoint() override = default;

   void set_initial_point(const std::vector<double>& initial_point) override;
   double compute_initial_value(double value, const Range& bounds) override;
   void set_constraints(const Problem& problem, Iterate& iterate);
   void initialize(Statistics& statistics, const Problem& problem, Iterate& first_iterate) override;
   void create_current_subproblem(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) override;
   void build_objective_model(const Problem& problem, Iterate& current_iterate, double objective_multiplier) override;
   void add_variable(size_t i, double current_value, const Range& bounds, double objective_term, size_t j, double jacobian_term) override;
   Direction solve(Statistics& statistics, const Problem& problem, Iterate& current_iterate) override;
   Direction compute_second_order_correction(const Problem& problem, Iterate& trial_iterate) override;
   [[nodiscard]] PredictedReductionModel generate_predicted_reduction_model(const Problem& problem, const Direction& direction) const override;
   void compute_progress_measures(const Problem& problem, Iterate& iterate) override;
   void register_accepted_iterate(Iterate& iterate) override;
   [[nodiscard]] int get_hessian_evaluation_count() const override;

private:
   // barrier parameter
   double barrier_parameter;
   const double tolerance;
   const std::unique_ptr<HessianModel<typename LinearSolverType::matrix_type>> hessian_model;
   typename LinearSolverType::matrix_type kkt_matrix;
   const std::unique_ptr<LinearSolver<typename LinearSolverType::matrix_type>> linear_solver;
   const InteriorPointParameters parameters;

   // data structures
   std::vector<size_t> lower_bounded_variables{}; // indices of the lower-bounded variables
   std::vector<size_t> upper_bounded_variables{}; // indices of the upper-bounded variables

   double inertia_hessian{0.};
   double inertia_hessian_last{0.};
   double inertia_constraints{0.};
   double default_multiplier;
   size_t iteration{0};
   size_t number_factorizations{0};

   // preallocated vectors
   std::vector<double> solution_IPM;
   std::vector<double> barrier_constraints;
   std::vector<double> rhs;
   std::vector<double> lower_delta_z;
   std::vector<double> upper_delta_z;

   void update_barrier_parameter(const Iterate& current_iterate);
   void set_variables_bounds(const Problem& problem, const Iterate& current_iterate, double trust_region_radius) override;
   void factorize(const Problem& problem, typename LinearSolverType::matrix_type& current_kkt_matrix);
   double compute_barrier_directional_derivative(const std::vector<double>& solution);
   double evaluate_barrier_function(const Problem& problem, Iterate& iterate);
   double primal_fraction_to_boundary(const Iterate& current_iterate, const std::vector<double>& ipm_solution, double tau);
   double dual_fraction_to_boundary(const Iterate& current_iterate, double tau);
   void assemble_kkt_matrix(Iterate& current_iterate);
   void modify_inertia(const Problem& problem, typename LinearSolverType::matrix_type& kkt_matrix, size_t size_first_block, size_t size_second_block);
   void generate_kkt_rhs(const Iterate& current_iterate);
   void compute_lower_bound_dual_direction(const Iterate& current_iterate, const std::vector<double>& solution);
   void compute_upper_bound_dual_direction(const Iterate& current_iterate, const std::vector<double>& solution);
   void generate_direction(const Problem& problem, const Iterate& current_iterate, std::vector<double>& solution_IPM);
   [[nodiscard]] double compute_KKT_error_scaling(const Iterate& current_iterate) const;
   [[nodiscard]] double compute_central_complementarity_error(const Iterate& iterate) const;
   void print_soc_iteration(const Direction& direction_soc) const;
};

template<typename LinearSolverType>
inline InteriorPoint<LinearSolverType>::InteriorPoint(const Problem& problem, size_t max_number_variables, size_t number_constraints,
      const std::string& hessian_evaluation_method, double initial_barrier_parameter, double default_multiplier, double tolerance,
      bool use_trust_region) :
      // add the slacks to the variables
      Subproblem(problem.number_variables + problem.inequality_constraints.size(), max_number_variables + problem.inequality_constraints.size(),
            number_constraints, SOC_UPON_REJECTION),
      barrier_parameter(initial_barrier_parameter), tolerance(tolerance),
      // if no trust region is used, the problem should be convexified. However, the inertia of the augmented matrix will be corrected later
      hessian_model(HessianModelFactory<typename LinearSolverType::matrix_type>::create(hessian_evaluation_method,
            this->max_number_variables, problem.hessian_maximum_number_nonzeros, false)),
      kkt_matrix(this->max_number_variables + number_constraints,
            problem.hessian_maximum_number_nonzeros + this->max_number_variables /* regularization */ +
            2 * this->max_number_variables /* diagonal barrier terms */ + this->max_number_variables * number_constraints /* Jacobian */),
      linear_solver(LinearSolverFactory<LinearSolverType>::create(this->max_number_variables + number_constraints)),
      parameters({0.99, 1e10, 100., 0.2, 1.5, 10.}), default_multiplier(default_multiplier),
      //current_primal_iterate(this->max_number_variables),
      //current_lower_bound_multipliers(this->max_number_variables), current_upper_bound_multipliers(this->max_number_variables),
      solution_IPM(this->max_number_variables + number_constraints),
      barrier_constraints(number_constraints), rhs(this->max_number_variables + number_constraints),
      lower_delta_z(this->max_number_variables), upper_delta_z(this->max_number_variables) {
   // register the original variables bounds
   copy_from(this->variables_bounds, problem.variables_bounds, problem.number_variables);

   // constraints are transformed into "c(x) = 0"
   for (size_t j = 0; j < problem.number_constraints; j++) {
      this->constraints_bounds[j] = {0., 0.};
   }

   // identify the bounded variables
   for (size_t i = 0; i < problem.number_variables; i++) {
      if (use_trust_region || (problem.variable_status[i] == BOUNDED_LOWER || problem.variable_status[i] == BOUNDED_BOTH_SIDES)) {
         this->lower_bounded_variables.push_back(i);
      }
      if (use_trust_region || (problem.variable_status[i] == BOUNDED_UPPER || problem.variable_status[i] == BOUNDED_BOTH_SIDES)) {
         this->upper_bounded_variables.push_back(i);
      }
   }
   // identify the inequality constraint slacks
   DEBUG << problem.inequality_constraints.size() << " slacks\n";
   for (const auto[j, i]: problem.inequality_constraints) {
      const size_t slack_index = problem.number_variables + i;
      if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
         this->lower_bounded_variables.push_back(slack_index);
      }
      if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
         this->upper_bounded_variables.push_back(slack_index);
      }
      // store the bounds of the slacks
      this->variables_bounds[slack_index] = problem.constraint_bounds[j];
   }
}

template<typename LinearSolverType>
inline void InteriorPoint<LinearSolverType>::set_initial_point(const std::vector<double>& /*initial_point*/) {
   // do nothing
}

template<typename LinearSolverType>
inline double InteriorPoint<LinearSolverType>::compute_initial_value(double value, const Range& bounds) {
   // the initial value of a variable must be strictly within bounds
   return Subproblem::push_variable_to_interior(value, bounds);
}

template<typename LinearSolverType>
inline void InteriorPoint<LinearSolverType>::set_constraints(const Problem& problem, Iterate& iterate) {
   iterate.evaluate_constraints(problem);
   // transform the constraints into "= 0" equalities
   for (const auto& element: problem.equality_constraints) {
      const size_t j = element.first;
      this->barrier_constraints[j] = iterate.constraints[j] - problem.constraint_bounds[j].lb;
   }
   for (const auto[j, i]: problem.inequality_constraints) {
      this->barrier_constraints[j] = iterate.constraints[j] - iterate.x[problem.number_variables + i];
   }
}

template<typename LinearSolverType>
inline void InteriorPoint<LinearSolverType>::initialize(Statistics& statistics, const Problem& problem, Iterate& first_iterate) {
   statistics.add_column("barrier param.", Statistics::double_width, 8);

   // resize to the new size (primals + slacks)
   first_iterate.change_number_variables(this->number_variables);

   // make the initial point strictly feasible wrt the bounds
   for (size_t i = 0; i < problem.number_variables; i++) {
      first_iterate.x[i] = this->compute_initial_value(first_iterate.x[i], problem.variables_bounds[i]);
   }

   // initialize the slacks and add contribution to the constraint Jacobian
   first_iterate.evaluate_constraints(problem);
   first_iterate.evaluate_constraint_jacobian(problem);
   for (const auto[j, i]: problem.inequality_constraints) {
      const double slack_value = this->compute_initial_value(first_iterate.constraints[j], problem.constraint_bounds[j]);
      first_iterate.x[problem.number_variables + i] = slack_value;
      first_iterate.constraint_jacobian[j].insert(problem.number_variables + i, -1.);
   }

   // compute least-square multipliers
   if (problem.is_constrained()) {
      Subproblem::compute_least_square_multipliers(problem, this->kkt_matrix, this->rhs, *this->linear_solver, first_iterate,
            first_iterate.multipliers.constraints);
   }

   // set the bound multipliers
   for (size_t i: this->lower_bounded_variables) {
      first_iterate.multipliers.lower_bounds[i] = this->default_multiplier;
   }
   for (size_t i: this->upper_bounded_variables) {
      first_iterate.multipliers.upper_bounds[i] = -this->default_multiplier;
   }

   // compute the optimality and feasibility measures of the initial point
   this->set_constraints(problem, first_iterate);
   this->compute_progress_measures(problem, first_iterate);
}

template<typename LinearSolverType>
inline void InteriorPoint<LinearSolverType>::create_current_subproblem(const Problem& problem, Iterate& current_iterate, double objective_multiplier,
      double trust_region_radius) {
   // update the barrier parameter if the current iterate solves the subproblem
   this->update_barrier_parameter(current_iterate);

   // constraints
   this->set_constraints(problem, current_iterate);
   copy_from(this->constraints_multipliers, current_iterate.multipliers.constraints);

   // constraint Jacobian
   problem.evaluate_constraint_jacobian(current_iterate.x, this->constraint_jacobian);
   // add the slack variables
   for (const auto[j, i]: problem.inequality_constraints) {
      this->constraint_jacobian[j].insert(problem.number_variables + i, -1.);
   }

   // build a model of the objective scaled by the objective multiplier
   this->build_objective_model(problem, current_iterate, objective_multiplier);

   // variables and bounds
   this->set_variables_bounds(problem, current_iterate, trust_region_radius);
}

template<typename LinearSolverType>
inline void InteriorPoint<LinearSolverType>::build_objective_model(const Problem& problem, Iterate& current_iterate, double objective_multiplier) {
   // evaluate the Hessian
   this->hessian_model->evaluate(problem, current_iterate.x, objective_multiplier, this->constraints_multipliers);

   // objective gradient
   this->set_scaled_objective_gradient(problem, current_iterate, objective_multiplier);
   for (size_t i: this->lower_bounded_variables) {
      this->objective_gradient.insert(i, -this->barrier_parameter / (current_iterate.x[i] - this->variables_bounds[i].lb));
   }
   for (size_t i: this->upper_bounded_variables) {
      this->objective_gradient.insert(i, -this->barrier_parameter / (current_iterate.x[i] - this->variables_bounds[i].ub));
   }
}

template<typename LinearSolverType>
inline void InteriorPoint<LinearSolverType>::add_variable(size_t i, double current_value, const Range& bounds, double objective_term, size_t j,
      double jacobian_term) {
   Subproblem::add_variable(i, current_value, bounds, objective_term, j, jacobian_term);
   if (-std::numeric_limits<double>::infinity() < bounds.lb) {
      this->lower_bounded_variables.push_back(i);
   }
   if (bounds.ub < std::numeric_limits<double>::infinity()) {
      this->upper_bounded_variables.push_back(i);
   }
}

template<typename LinearSolverType>
inline Direction InteriorPoint<LinearSolverType>::solve(Statistics& statistics, const Problem& problem, Iterate& current_iterate) {
   this->iteration++;
   // assemble and factorize the KKT matrix
   this->assemble_kkt_matrix(current_iterate);
   this->factorize(problem, this->kkt_matrix);
   // inertia correction
   this->modify_inertia(problem, this->kkt_matrix, this->number_variables, problem.number_constraints);
   DEBUG << "KKT matrix:\n" << this->kkt_matrix << "\n";
   auto[number_pos, number_neg, number_zero] = this->linear_solver->get_inertia();
   assert(number_pos == this->number_variables && number_neg == problem.number_constraints && number_zero == 0);

   // right-hand side
   this->generate_kkt_rhs(current_iterate);

   // compute the solution (Δx, -Δλ)
   this->linear_solver->solve(this->kkt_matrix, this->rhs, this->solution_IPM);
   this->number_subproblems_solved++;

   // generate IPM direction
   this->generate_direction(problem, current_iterate, this->solution_IPM);
   statistics.add_statistic("barrier param.", this->barrier_parameter);
   return this->direction;
}

//   catch (const UnstableInertiaCorrection& e) {
//      // unstable factorization during optimality phase
//      throw "InteriorPoint: inertia correction failed";
//   }

template<typename LinearSolverType>
inline Direction InteriorPoint<LinearSolverType>::compute_second_order_correction(const Problem& problem, Iterate& trial_iterate) {
   // modify the RHS by adding the values of the constraints
   for (const auto& element: problem.equality_constraints) {
      size_t j = element.first;
      this->rhs[this->number_variables + j] -= trial_iterate.constraints[j] - problem.constraint_bounds[j].lb;
   }
   for (const auto[j, i]: problem.inequality_constraints) {
      this->rhs[this->number_variables + j] -= trial_iterate.constraints[j] - trial_iterate.x[problem.number_variables + i];
   }

   // compute the solution (Δx, -Δλ)
   this->linear_solver->solve(kkt_matrix, this->rhs, this->solution_IPM);
   this->number_subproblems_solved++;

   // generate IPM direction
   this->generate_direction(problem, trial_iterate, this->solution_IPM);
   this->print_soc_iteration(this->direction);
   return this->direction;
}

template<typename LinearSolverType>
void InteriorPoint<LinearSolverType>::print_soc_iteration(const Direction& direction_soc) const {
   DEBUG << "Entered SOC computation\n";
   DEBUG << "KKT matrix:\n" << this->kkt_matrix << "\n";
   DEBUG << "SOC RHS: ";
   print_vector(DEBUG, this->rhs);
   DEBUG << "\n";
   DEBUG << "SOC direction:\n" << direction_soc << "\n";
}

template<typename LinearSolverType>
inline PredictedReductionModel InteriorPoint<LinearSolverType>::generate_predicted_reduction_model(const Problem& /*problem*/,
      const Direction& direction) const {
   return PredictedReductionModel(-direction.objective, [&]() {
      return [=](double step_length) {
         return -step_length * direction.objective;
      };
   });
}

template<typename LinearSolverType>
inline void InteriorPoint<LinearSolverType>::compute_progress_measures(const Problem& problem, Iterate& iterate) {
   auto residual_function = [&](size_t j) {
      if (problem.constraint_status[j] == EQUAL_BOUNDS) {
         return iterate.constraints[j] - problem.constraint_bounds[j].lb;
      }
      else {
         const size_t i = problem.inequality_constraints.at(j);
         return iterate.constraints[j] - iterate.x[problem.number_variables + i];
      }
   };
   const double constraint_violation = norm_1(residual_function, problem.number_constraints);

   // compute barrier objective
   const double objective = this->evaluate_barrier_function(problem, iterate);
   iterate.progress = {constraint_violation, objective};
}

template<typename LinearSolverType>
inline void InteriorPoint<LinearSolverType>::register_accepted_iterate(Iterate& iterate) {
   // rescale the bound multipliers (Eq (16) in Ipopt paper)
   for (size_t i: this->lower_bounded_variables) {
      const double coefficient = this->barrier_parameter / (iterate.x[i] - this->variables_bounds[i].lb);
      const double lb = coefficient / this->parameters.k_sigma;
      const double ub = coefficient * this->parameters.k_sigma;
      assert(lb <= ub && "IPM bound multiplier reset: the bounds are in the wrong order");
      iterate.multipliers.lower_bounds[i] = std::max(std::min(iterate.multipliers.lower_bounds[i], ub), lb);
   }
   for (size_t i: this->upper_bounded_variables) {
      const double coefficient = this->barrier_parameter / (iterate.x[i] - this->variables_bounds[i].ub);
      const double lb = coefficient * this->parameters.k_sigma;
      const double ub = coefficient / this->parameters.k_sigma;
      assert(lb <= ub && "IPM bound multiplier reset: the bounds are in the wrong order");
      iterate.multipliers.upper_bounds[i] = std::max(std::min(iterate.multipliers.upper_bounds[i], ub), lb);
   }
}

template<typename LinearSolverType>
inline int InteriorPoint<LinearSolverType>::get_hessian_evaluation_count() const {
   return this->hessian_model->evaluation_count;
}

template<typename LinearSolverType>
inline void InteriorPoint<LinearSolverType>::update_barrier_parameter(const Iterate& current_iterate) {
   const double tolerance_fraction = this->tolerance / 10.;
   // scaled error terms
   const double sd = this->compute_KKT_error_scaling(current_iterate);
   const double KKTerror = current_iterate.errors.KKT / sd;
   const double central_complementarity_error = this->compute_central_complementarity_error(current_iterate);
   const double error = std::max({KKTerror, current_iterate.errors.constraints, central_complementarity_error});

   // update of the barrier problem
   while (error <= this->parameters.k_epsilon * this->barrier_parameter && tolerance_fraction < this->barrier_parameter) {
      this->barrier_parameter = std::max(tolerance_fraction, std::min(this->parameters.k_mu * this->barrier_parameter,
            std::pow(this->barrier_parameter, this->parameters.theta_mu)));
      DEBUG << "IPM: mu updated to " << this->barrier_parameter << " and filter reset\n";
      // signal the redefinition of the problem to the globalization strategy
      this->subproblem_definition_changed = true;
   }
   DEBUG << "mu is " << this->barrier_parameter << "\n";
}

template<typename LinearSolverType>
inline void InteriorPoint<LinearSolverType>::set_variables_bounds(const Problem& problem, const Iterate& current_iterate,
      double trust_region_radius) {
   // here, we work with the original bounds
   // very important: apply the trust region only on the original variables (not the slacks)
   for (size_t i = 0; i < problem.number_variables; i++) {
      double lb = std::max(current_iterate.x[i] - trust_region_radius, problem.variables_bounds[i].lb);
      double ub = std::min(current_iterate.x[i] + trust_region_radius, problem.variables_bounds[i].ub);
      this->variables_bounds[i] = {lb, ub};
   }
}

template<typename LinearSolverType>
inline void InteriorPoint<LinearSolverType>::factorize(const Problem& problem, typename LinearSolverType::matrix_type& current_kkt_matrix) {
   // compute the symbolic factorization only when:
   // the problem has a non-constant augmented system (ie is not an LP or a QP) or it is the first factorization
   if (this->number_factorizations == 0 || !problem.fixed_hessian_sparsity || problem.type == NONLINEAR) {
      this->linear_solver->do_symbolic_factorization(current_kkt_matrix);
   }
   this->linear_solver->do_numerical_factorization(current_kkt_matrix);
   this->number_factorizations++;
}

template<typename LinearSolverType>
inline double InteriorPoint<LinearSolverType>::compute_barrier_directional_derivative(const std::vector<double>& solution) {
   return dot(solution, this->objective_gradient);
}

template<typename LinearSolverType>
inline double InteriorPoint<LinearSolverType>::evaluate_barrier_function(const Problem& problem, Iterate& iterate) {
   double objective = 0.;
   // bound constraints
   for (size_t i: this->lower_bounded_variables) {
      objective -= std::log(iterate.x[i] - this->variables_bounds[i].lb);
   }
   for (size_t i: this->upper_bounded_variables) {
      objective -= std::log(this->variables_bounds[i].ub - iterate.x[i]);
   }
   objective *= this->barrier_parameter;
   // original objective
   iterate.evaluate_objective(problem);
   objective += iterate.objective;
   return objective;
}

template<typename LinearSolverType>
inline double InteriorPoint<LinearSolverType>::primal_fraction_to_boundary(const Iterate& current_iterate, const std::vector<double>& ipm_solution,
      double tau) {
   double primal_length = 1.;
   for (size_t i: this->lower_bounded_variables) {
      if (ipm_solution[i] < 0.) {
         double trial_alpha_xi = -tau * (current_iterate.x[i] - this->variables_bounds[i].lb) / ipm_solution[i];
         primal_length = std::min(primal_length, trial_alpha_xi);
      }
   }
   for (size_t i: this->upper_bounded_variables) {
      if (0. < ipm_solution[i]) {
         double trial_alpha_xi = -tau * (current_iterate.x[i] - this->variables_bounds[i].ub) / ipm_solution[i];
         primal_length = std::min(primal_length, trial_alpha_xi);
      }
   }
   return primal_length;
}

template<typename LinearSolverType>
inline double InteriorPoint<LinearSolverType>::dual_fraction_to_boundary(const Iterate& current_iterate, double tau) {
   double dual_length = 1.;
   for (size_t i = 0; i < this->number_variables; i++) {
      if (this->lower_delta_z[i] < 0.) {
         double trial_alpha_zj = -tau * current_iterate.multipliers.lower_bounds[i] / this->lower_delta_z[i];
         dual_length = std::min(dual_length, trial_alpha_zj);
      }
      if (0. < this->upper_delta_z[i]) {
         double trial_alpha_zj = -tau * current_iterate.multipliers.upper_bounds[i] / this->upper_delta_z[i];
         dual_length = std::min(dual_length, trial_alpha_zj);
      }
   }
   return dual_length;
}

template<typename LinearSolverType>
inline void InteriorPoint<LinearSolverType>::assemble_kkt_matrix(Iterate& current_iterate) {
   this->kkt_matrix.reset();
   // copy the Lagrangian Hessian
   // assume that the Hessian is sorted
   size_t current_column = 0;
   this->hessian_model->hessian.for_each([&](int i, int j, double entry) {
      if (j != static_cast<int>(current_column)) {
         for (size_t column = current_column; column < static_cast<size_t>(j); column++) {
            this->kkt_matrix.finalize(column);
            current_column++;
         }
      }
      this->kkt_matrix.insert(entry, i, j);
   });

   // diagonal terms: bounds of primals and slacks
   for (size_t i: this->lower_bounded_variables) {
      this->kkt_matrix.insert(current_iterate.multipliers.lower_bounds[i] / (current_iterate.x[i] - this->variables_bounds[i].lb), i, i);
   }
   for (size_t i: this->upper_bounded_variables) {
      this->kkt_matrix.insert(current_iterate.multipliers.upper_bounds[i] / (current_iterate.x[i] - this->variables_bounds[i].ub), i, i);
   }

   // Jacobian of general constraints
   for (size_t j = 0; j < this->number_constraints; j++) {
      this->constraint_jacobian[j].for_each([&](size_t i, double derivative) {
         this->kkt_matrix.insert(derivative, i, this->number_variables + j);
      });
      this->kkt_matrix.finalize(j);
   }
}

template<typename LinearSolverType>
inline void InteriorPoint<LinearSolverType>::modify_inertia(const Problem& problem, typename LinearSolverType::matrix_type& kkt_matrix, size_t
size_first_block, size_t size_second_block) {
   this->inertia_hessian = 0.;
   this->inertia_constraints = 0.;
   DEBUG << "Testing factorization with inertia term " << this->inertia_hessian << "\n";

   bool good_inertia = false;
   if (!this->linear_solver->matrix_is_singular() && this->linear_solver->number_negative_eigenvalues() == size_second_block) {
      DEBUG << "Inertia is good\n";
      good_inertia = true;
   }
   else {
      DEBUG << "Inertia is not good\n";
      // inertia term for constraints
      if (this->linear_solver->matrix_is_singular()) {
         DEBUG << "Matrix is singular\n";
         this->inertia_constraints = 1e-8 * std::pow(this->barrier_parameter, 0.25);
      }
      else {
         this->inertia_constraints = 0.;
      }
      // inertia term for Hessian
      if (this->inertia_hessian_last == 0.) {
         this->inertia_hessian = 1e-4;
      }
      else {
         this->inertia_hessian = std::max(1e-20, this->inertia_hessian_last / 3.);
      }
   }

   size_t current_matrix_size = kkt_matrix.number_nonzeros;
   if (!good_inertia) {
      for (size_t i = 0; i < size_first_block; i++) {
         kkt_matrix.insert(this->inertia_hessian, i, i);
      }
      for (size_t j = size_first_block; j < size_first_block + size_second_block; j++) {
         kkt_matrix.insert(-this->inertia_constraints, j, j);
      }
   }

   while (!good_inertia) {
      DEBUG << "Testing factorization with inertia term " << this->inertia_hessian << "\n";
      DEBUG << kkt_matrix << "\n";
      this->factorize(problem, kkt_matrix);

      if (!this->linear_solver->matrix_is_singular() && this->linear_solver->number_negative_eigenvalues() == size_second_block) {
         good_inertia = true;
         DEBUG << "Factorization was a success\n";
         this->inertia_hessian_last = this->inertia_hessian;
      }
      else {
         if (this->inertia_hessian_last == 0.) {
            this->inertia_hessian *= 100.;
         }
         else {
            this->inertia_hessian *= 8.;
         }

         if (this->inertia_hessian <= 1e40) {
            for (size_t i = 0; i < size_first_block; i++) {
               kkt_matrix.matrix[current_matrix_size + i] = this->inertia_hessian;
            }
            for (size_t j = size_first_block; j < size_first_block + size_second_block; j++) {
               kkt_matrix.matrix[current_matrix_size + j] = -this->inertia_constraints;
            }
         }
         else {
            throw UnstableInertiaCorrection();
         }
      }
   }
}

template<typename LinearSolverType>
inline void InteriorPoint<LinearSolverType>::generate_kkt_rhs(const Iterate& current_iterate) {
   // generate the right-hand side
   clear(this->rhs);

   // barrier objective gradient
   this->objective_gradient.for_each([&](size_t i, double derivative) {
      this->rhs[i] = -derivative;
   });

   // constraint: evaluations and gradients
   for (size_t j = 0; j < current_iterate.constraints.size(); j++) {
      // Lagrangian
      if (this->constraints_multipliers[j] != 0.) {
         this->constraint_jacobian[j].for_each([&](size_t i, double derivative) {
            this->rhs[i] += this->constraints_multipliers[j] * derivative;
         });
      }
      // constraints
      this->rhs[this->number_variables + j] = -this->barrier_constraints[j];
   }
   DEBUG << "RHS: ";
   print_vector(DEBUG, this->rhs);
   DEBUG << "\n";
}

template<typename LinearSolverType>
inline void InteriorPoint<LinearSolverType>::compute_lower_bound_dual_direction(const Iterate& current_iterate, const std::vector<double>&
solution) {
   clear(this->lower_delta_z);
   for (size_t i: this->lower_bounded_variables) {
      const double distance_to_bound = current_iterate.x[i] - this->variables_bounds[i].lb;
      this->lower_delta_z[i] = (this->barrier_parameter - solution[i] * current_iterate.multipliers.lower_bounds[i]) / distance_to_bound -
                               current_iterate.multipliers.lower_bounds[i];
   }
}

template<typename LinearSolverType>
inline void InteriorPoint<LinearSolverType>::compute_upper_bound_dual_direction(const Iterate& current_iterate, const std::vector<double>&
solution) {
   clear(this->upper_delta_z);
   for (size_t i: this->upper_bounded_variables) {
      const double distance_to_bound = current_iterate.x[i] - this->variables_bounds[i].ub;
      this->upper_delta_z[i] = (this->barrier_parameter - solution[i] * current_iterate.multipliers.upper_bounds[i]) / distance_to_bound -
                               current_iterate.multipliers.upper_bounds[i];
   }
}

template<typename LinearSolverType>
inline void InteriorPoint<LinearSolverType>::generate_direction(const Problem& problem, const Iterate& current_iterate,
      std::vector<double>& solution_IPM) {
   // retrieve +Δλ (Nocedal p590)
   for (size_t j = this->number_variables; j < solution_IPM.size(); j++) {
      solution_IPM[j] = -solution_IPM[j];
   }

   // "fraction to boundary" rule for primal variables and constraints multipliers
   const double tau = std::max(this->parameters.tau_min, 1. - this->barrier_parameter);
   const double primal_step_length = this->primal_fraction_to_boundary(current_iterate, solution_IPM, tau);
   for (size_t i = 0; i < this->number_variables; i++) {
      this->direction.x[i] = primal_step_length * solution_IPM[i];
   }
   for (size_t j = 0; j < problem.number_constraints; j++) {
      this->direction.multipliers.constraints[j] = primal_step_length * solution_IPM[this->number_variables + j];
   }

   // compute bound multiplier direction Δz
   this->compute_lower_bound_dual_direction(current_iterate, solution_IPM);
   this->compute_upper_bound_dual_direction(current_iterate, solution_IPM);

   // "fraction to boundary" rule for bound multipliers
   const double dual_step_length = this->dual_fraction_to_boundary(current_iterate, tau);
   for (size_t i = 0; i < this->number_variables; i++) {
      this->direction.multipliers.lower_bounds[i] = current_iterate.multipliers.lower_bounds[i] + dual_step_length * this->lower_delta_z[i];
      this->direction.multipliers.upper_bounds[i] = current_iterate.multipliers.upper_bounds[i] + dual_step_length * this->upper_delta_z[i];
   }

   this->direction.status = OPTIMAL;
   this->direction.norm = norm_inf(direction.x, 0, this->number_variables);
   // evaluate the barrier objective
   this->direction.objective = this->compute_barrier_directional_derivative(direction.x);

   DEBUG << "IPM solution:\n";
   DEBUG << "Δx: ";
   print_vector(DEBUG, solution_IPM, '\n', 0, problem.number_variables);
   DEBUG << "Δs: ";
   print_vector(DEBUG, solution_IPM, '\n', problem.number_variables, problem.inequality_constraints.size());
   DEBUG << "Δλ: ";
   print_vector(DEBUG, solution_IPM, '\n', this->number_variables, problem.number_constraints);
   DEBUG << "Δz_L: ";
   print_vector(DEBUG, this->lower_delta_z);
   DEBUG << "Δz_U: ";
   print_vector(DEBUG, this->upper_delta_z);
   DEBUG << "primal length = " << primal_step_length << "\n";
   DEBUG << "dual length = " << dual_step_length << "\n\n";
}

template<typename LinearSolverType>
inline double InteriorPoint<LinearSolverType>::compute_KKT_error_scaling(const Iterate& current_iterate) const {
   // KKT error
   const double norm_1_constraint_multipliers = norm_1(current_iterate.multipliers.constraints);
   const double norm_1_bound_multipliers = norm_1(current_iterate.multipliers.lower_bounds) + norm_1(current_iterate.multipliers.upper_bounds);
   const double norm_1_multipliers = norm_1_constraint_multipliers + norm_1_bound_multipliers;
   const size_t total_size = this->number_variables + current_iterate.multipliers.constraints.size();
   const double sd = std::max(this->parameters.smax, norm_1_multipliers / static_cast<double>(total_size)) / this->parameters.smax;
   return sd;
}

template<typename LinearSolverType>
inline double InteriorPoint<LinearSolverType>::compute_central_complementarity_error(const Iterate& iterate) const {
   // variable bound constraints
   const auto residual_function = [&](size_t i) {
      double result = 0.;
      if (-std::numeric_limits<double>::infinity() < this->variables_bounds[i].lb) {
         result += iterate.multipliers.lower_bounds[i] * (iterate.x[i] - this->variables_bounds[i].lb) - this->barrier_parameter;
      }
      if (this->variables_bounds[i].ub < std::numeric_limits<double>::infinity()) {
         result += iterate.multipliers.upper_bounds[i] * (iterate.x[i] - this->variables_bounds[i].ub) - this->barrier_parameter;
      }
      return result;
   };

   // scaling
   const double bound_multipliers_norm = norm_1(iterate.multipliers.lower_bounds) + norm_1(iterate.multipliers.upper_bounds);
   const double sc = std::max(this->parameters.smax, bound_multipliers_norm / static_cast<double>(this->number_variables)) / this->parameters.smax;
   return norm_1(residual_function, this->number_variables) / sc;
}

#endif // IPM_H
