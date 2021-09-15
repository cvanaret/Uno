#include "Preprocessing.hpp"
#include "solvers/QP/BQPDSolver.hpp"

void Preprocessing::apply(const Problem& problem, Iterate& first_iterate) {
   /* linear constraints */
   INFO << "Preprocessing phase: the problem has " << problem.linear_constraints.size() << " linear constraints\n";
   if (!problem.linear_constraints.empty()) {
      first_iterate.compute_constraints(problem);

      int infeasible_linear_constraints = 0;
      for (const auto element: problem.linear_constraints) {
         size_t j = element.first;
         if (first_iterate.constraints[j] < problem.constraint_bounds[j].lb || problem.constraint_bounds[j].ub < first_iterate.constraints[j]) {
            infeasible_linear_constraints++;
         }
      }
      INFO << "There are " << infeasible_linear_constraints << " infeasible linear constraints at the initial point\n";

      if (0 < infeasible_linear_constraints) {
         INFO << "Current point: ";
         print_vector(INFO, first_iterate.x);
         size_t number_constraints = problem.linear_constraints.size();
         BQPDSolver solver(problem.number_variables, number_constraints, problem.number_variables, true);

         CSCSymmetricMatrix hessian = CSCSymmetricMatrix::identity(problem.number_variables);
         SparseVector<double> linear_objective(0); // empty
         // constraints Jacobian
         std::vector<SparseVector<double>> constraints_jacobian(number_constraints);
         for (size_t j = 0; j < problem.number_constraints; j++) {
            constraints_jacobian[j].reserve(problem.number_variables);
         }
         for (const auto[j, linear_constraint_index]: problem.linear_constraints) {
            problem.evaluate_constraint_gradient(first_iterate.x, j, constraints_jacobian[linear_constraint_index]);
         }
         // variables bounds
         std::vector<Range> variables_bounds(problem.number_variables);
         for (size_t i = 0; i < problem.number_variables; i++) {
            variables_bounds[i] = {problem.variables_bounds[i].lb - first_iterate.x[i], problem.variables_bounds[i].ub - first_iterate.x[i]};
         }
         // constraints bounds
         std::vector<Range> constraints_bounds(number_constraints);
         for (const auto[j, linear_constraint_index]: problem.linear_constraints) {
            constraints_bounds[linear_constraint_index] =
                  {problem.constraint_bounds[j].lb - first_iterate.constraints[j], problem.constraint_bounds[j].ub - first_iterate.constraints[j]};
         }
         std::vector<double> d0(problem.number_variables);
         Direction direction = solver.solve_QP(variables_bounds, constraints_bounds, linear_objective, constraints_jacobian, hessian, d0);
         if (direction.status == INFEASIBLE) {
            throw std::runtime_error("Linear constraints cannot be satisfied");
         }

         add_vectors(first_iterate.x, direction.x, 1., first_iterate.x);
         // copy bound multipliers
         first_iterate.multipliers.lower_bounds = direction.multipliers.lower_bounds;
         first_iterate.multipliers.upper_bounds = direction.multipliers.upper_bounds;
         // copy constraint multipliers
         for (const auto[j, linear_constraint_index]: problem.linear_constraints) {
            first_iterate.multipliers.constraints[j] = direction.multipliers.constraints[linear_constraint_index];
         }
         INFO << "Linear feasible initial point: ";
         print_vector(INFO, first_iterate.x);
         INFO << "\n";
      }
   }
}