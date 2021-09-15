#ifndef SLP_H
#define SLP_H

#include "Subproblem.hpp"
#include "solvers/QP/LPSolver.hpp"
#include "solvers/QP/LPSolverFactory.hpp"

template<typename LPSolverType>
class SLP : public Subproblem {
public:
   SLP(const Problem& problem, size_t number_variables, size_t number_constraints);

   void create_current_subproblem(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) override;
   void build_objective_model(const Problem& problem, Iterate& current_iterate, double objective_multiplier) override;
   void set_initial_point(const std::vector<double>& point) override;
   Direction solve(Statistics& statistics, const Problem& problem, Iterate& current_iterate) override;
   [[nodiscard]] PredictedReductionModel generate_predicted_reduction_model(const Problem& problem, const Direction& direction) const override;
   [[nodiscard]] int get_hessian_evaluation_count() const override;

private:
   /* use pointers to allow polymorphism */
   const std::unique_ptr <LPSolverType> solver; /*!< Solver that solves the subproblem */
   std::vector<double> initial_point;
};

template<typename LPSolverType>
inline SLP<LPSolverType>::SLP(const Problem& problem, size_t number_variables, size_t number_constraints) : Subproblem(number_variables,
      number_constraints),
      solver(LPSolverFactory<LPSolverType>::create(number_variables, number_constraints)),
      initial_point(number_variables) {
   // register the original constraints bounds
   for (size_t j = 0; j < problem.number_constraints; j++) {
      this->constraints_bounds[j] = problem.constraint_bounds[j];
   }
}

template<typename LPSolverType>
inline void SLP<LPSolverType>::create_current_subproblem(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) {
   copy_from(this->constraints_multipliers, current_iterate.multipliers.constraints);
   /* compute first- and second-order information */
   problem.evaluate_constraints(current_iterate.x, current_iterate.constraints);
   for (auto& row: this->constraints_jacobian) {
      row.clear();
   }
   problem.evaluate_constraints_jacobian(current_iterate.x, this->constraints_jacobian);

   this->build_objective_model(problem, current_iterate, objective_multiplier);

   /* bounds of the variables */
   this->set_variables_bounds(problem, current_iterate, trust_region_radius);

   /* bounds of the linearized constraints */
   this->set_constraints_bounds(problem, current_iterate.constraints);

   /* set the initial point */
   clear(this->initial_point);
}

template<typename LPSolverType>
inline void SLP<LPSolverType>::build_objective_model(const Problem& problem, Iterate& current_iterate, double objective_multiplier) {
   // objective gradient
   this->set_scaled_objective_gradient(problem, current_iterate, objective_multiplier);

   // initial point
   clear(this->initial_point);
}

template<typename LPSolverType>
inline void SLP<LPSolverType>::set_initial_point(const std::vector<double>& point) {
   copy_from(this->initial_point, point);
}

template<typename LPSolverType>
inline Direction SLP<LPSolverType>::solve(Statistics& /*statistics*/, const Problem& problem, Iterate& current_iterate) {
   /* solve the LP */
   Direction direction = this->solver->solve_LP(variables_bounds, constraints_bounds, this->objective_gradient, this->constraints_jacobian,
         this->initial_point);
   // compute dual displacements (SQP methods compute the new duals, not the displacements)
   for (size_t j = 0; j < problem.number_constraints; j++) {
      direction.multipliers.constraints[j] -= current_iterate.multipliers.constraints[j];
   }
   this->number_subproblems_solved++;
   return direction;
}

template<typename QPSolverType>
inline PredictedReductionModel SLP<QPSolverType>::generate_predicted_reduction_model(const Problem& /*problem*/, const Direction& direction) const {
   return PredictedReductionModel(-direction.objective, [&]() { // capture direction by reference
      // return a function of the step length that cheaply assembles the predicted reduction
      return [=](double step_length) { // capture the expensive quantities by value
         return -step_length * direction.objective;
      };
   });
}

template<typename LPSolverType>
inline int SLP<LPSolverType>::get_hessian_evaluation_count() const {
   // no second order evaluation is used
   return 0;
}

#endif // SLP_H
