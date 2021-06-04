#include "Preprocessing.hpp"
#include "BQPDSolver.hpp"

void Preprocessing::apply(Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
   /* linear constraints */
   INFO << "Preprocessing phase: the problem has " << problem.linear_constraints.size() << " linear constraints\n";
   if (!problem.linear_constraints.empty()) {
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
         INFO << "Current point: ";
         print_vector(INFO, x);
         int number_constraints = (int) problem.linear_constraints.size();
         BQPDSolver solver(problem.number_variables, number_constraints, problem.number_variables, true);

         int fortran_indexing = 1;
         CSCMatrix hessian = CSCMatrix::identity(problem.number_variables, fortran_indexing);
         SparseVector linear_objective; // empty
         std::vector<double> d0(problem.number_variables);
         // constraints Jacobian
         std::vector<SparseVector> constraints_jacobian(number_constraints);
         for (const auto[j, linear_constraint_index]: problem.linear_constraints) {
            problem.constraint_gradient(x, j, constraints_jacobian[linear_constraint_index]);
         }
         // variables bounds
         std::vector<Range> variables_bounds(problem.number_variables);
         for (size_t i = 0; i < problem.number_variables; i++) {
            variables_bounds[i] = {problem.variables_bounds[i].lb - x[i], problem.variables_bounds[i].ub - x[i]};
         }
         // constraints bounds
         std::vector<Range> constraints_bounds(number_constraints);
         for (const auto[j, linear_constraint_index]: problem.linear_constraints) {
            constraints_bounds[linear_constraint_index] =
                  {problem.constraint_bounds[j].lb - constraints[j], problem.constraint_bounds[j].ub - constraints[j]};
         }
         Direction direction = solver.solve_QP(variables_bounds, constraints_bounds, linear_objective, constraints_jacobian, hessian, d0);
         if (direction.status == INFEASIBLE) {
            throw std::runtime_error("Linear constraints cannot be satisfied");
         }

         std::vector<double> feasible_x = add_vectors(x, direction.x, 1.);
         x = feasible_x;
         // copy bound multipliers
         multipliers.lower_bounds = direction.multipliers.lower_bounds;
         multipliers.upper_bounds = direction.multipliers.upper_bounds;
         // copy constraint multipliers
         for (const auto[j, linear_constraint_index]: problem.linear_constraints) {
            multipliers.constraints[j] = direction.multipliers.constraints[linear_constraint_index];
         }
         INFO << "Linear feasible initial point: ";
         print_vector(INFO, x);
         INFO << "\n";
      }
   }
}