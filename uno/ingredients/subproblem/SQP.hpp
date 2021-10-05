#ifndef SQP_H
#define SQP_H

#include <cassert>
#include "Subproblem.hpp"
#include "HessianModel.hpp"
#include "solvers/QP/QPSolver.hpp"
#include "solvers/QP/QPSolverFactory.hpp"

template<typename QPSolverType>
class SQP : public Subproblem {
public:
   SQP(const Problem& problem, size_t max_number_variables, size_t number_constraints, const std::string& hessian_evaluation_method,
         bool use_trust_region);

   void create_current_subproblem(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) override;
   void build_objective_model(const Problem& problem, Iterate& current_iterate, double objective_multiplier) override;
   void set_initial_point(const std::vector<double>& point) override;
   Direction solve(Statistics& statistics, const Problem& problem, Iterate& current_iterate) override;
   [[nodiscard]] PredictedReductionModel generate_predicted_reduction_model(const Problem& problem, const Direction& direction) const override;
   [[nodiscard]] int get_hessian_evaluation_count() const override;

protected:
   // use pointers to allow polymorphism
   const std::unique_ptr <QPSolver<typename QPSolverType::matrix_type>> solver; /*!< Solver that solves the subproblem */
   const std::unique_ptr <HessianModel<typename QPSolverType::matrix_type>> hessian_model; /*!< Strategy to evaluate or approximate the
 * Hessian */
   std::vector<double> initial_point;
};

template<typename QPSolverType>
inline SQP<QPSolverType>::SQP(const Problem& problem, size_t max_number_variables, size_t number_constraints,
      const std::string& hessian_evaluation_method, bool use_trust_region) :
      Subproblem(problem.number_variables, max_number_variables, number_constraints, NO_SOC),
      // maximum number of Hessian nonzeros = number nonzeros + possible diagonal inertia correction
      solver(QPSolverFactory<QPSolverType>::create(max_number_variables, number_constraints,
            problem.hessian_maximum_number_nonzeros + max_number_variables, true)),
      // if no trust region is used, the problem should be convexified by controlling the inertia of the Hessian
      hessian_model(HessianModelFactory<typename QPSolverType::matrix_type>::create(hessian_evaluation_method, max_number_variables,
            problem.hessian_maximum_number_nonzeros + problem.number_variables, !use_trust_region)),
      initial_point(max_number_variables) {
}

template<typename QPSolverType>
inline void SQP<QPSolverType>::create_current_subproblem(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) {
   copy_from(this->constraints_multipliers, current_iterate.multipliers.constraints);

   // constraints
   current_iterate.evaluate_constraints(problem);

   // constraint Jacobian
   current_iterate.evaluate_constraint_jacobian(problem);
   this->constraint_jacobian = current_iterate.constraint_jacobian;

   // objective
   this->build_objective_model(problem, current_iterate, objective_multiplier);

   // bounds of the variables
   this->set_variables_bounds(problem, current_iterate, trust_region_radius);

   // bounds of the linearized constraints
   this->set_constraints_bounds(problem, current_iterate.constraints);

   // reset the initial point
   clear(this->initial_point);
}

template<typename QPSolverType>
inline void SQP<QPSolverType>::build_objective_model(const Problem& problem, Iterate& current_iterate, double objective_multiplier) {
   // Hessian
   this->hessian_model->evaluate(problem, current_iterate.x, objective_multiplier, this->constraints_multipliers);

   // objective gradient
   this->set_scaled_objective_gradient(problem, current_iterate, objective_multiplier);

   // initial point
   clear(this->initial_point);
}

template<typename QPSolverType>
inline void SQP<QPSolverType>::set_initial_point(const std::vector<double>& point) {
   copy_from(this->initial_point, point);
}

template<typename QPSolverType>
inline Direction SQP<QPSolverType>::solve(Statistics& /*statistics*/, const Problem& problem, Iterate& current_iterate) {
   /* compute QP direction */
   Direction direction = this->solver->solve_QP(this->variables_bounds, this->constraints_bounds, this->objective_gradient,
         this->constraint_jacobian, this->hessian_model->hessian, this->initial_point);
   this->number_subproblems_solved++;

   // compute dual displacements (SQP methods usually compute the new duals, not the displacements)
   for (size_t j = 0; j < problem.number_constraints; j++) {
      direction.multipliers.constraints[j] -= current_iterate.multipliers.constraints[j];
   }
   return direction;
}

template<typename QPSolverType>
inline PredictedReductionModel SQP<QPSolverType>::generate_predicted_reduction_model(const Problem& /*problem*/, const Direction& direction) const {
   return PredictedReductionModel(-direction.objective, [&]() { // capture this and direction by reference
      // precompute expensive quantities
      const double linear_term = dot(direction.x, this->objective_gradient);
      const double quadratic_term = this->hessian_model->hessian.quadratic_product(direction.x, direction.x) / 2.;
      // return a function of the step length that cheaply assembles the predicted reduction
      return [=](double step_length) { // capture the expensive quantities by value
         return -step_length * (linear_term + step_length * quadratic_term);
      };
   });
}

template<typename QPSolverType>
inline int SQP<QPSolverType>::get_hessian_evaluation_count() const {
   return this->hessian_model->evaluation_count;
}

#endif // SQP_H
