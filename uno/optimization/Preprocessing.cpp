#include "Preprocessing.hpp"
#include "solvers/QP/BQPDSolver.hpp"

void Preprocessing::enforce_linear_constraints(const Problem& problem, Iterate& first_iterate) {
   /* linear constraints */
   INFO << "Preprocessing phase: the problem has " << problem.linear_constraints.size() << " linear constraints\n";
   if (!problem.linear_constraints.empty()) {
      first_iterate.evaluate_constraints(problem);

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
         std::vector<SparseVector<double>> constraint_jacobian(number_constraints);
         for (size_t j = 0; j < problem.number_constraints; j++) {
            constraint_jacobian[j].reserve(problem.number_variables);
         }
         for (const auto[j, linear_constraint_index]: problem.linear_constraints) {
            problem.evaluate_constraint_gradient(first_iterate.x, j, constraint_jacobian[linear_constraint_index]);
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
         Direction direction = solver.solve_QP(variables_bounds, constraints_bounds, linear_objective, constraint_jacobian, hessian, d0);
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

// compute a least-square approximation of the multipliers by solving a linear system (uses existing linear system)
void Preprocessing::compute_least_square_multipliers(const Problem& problem, SymmetricMatrix& matrix, std::vector<double>& rhs, LinearSolver& solver,
      Iterate& current_iterate, std::vector<double>& multipliers, double multipliers_max_size) {
   const size_t number_variables = current_iterate.x.size();
   current_iterate.evaluate_objective_gradient(problem);
   current_iterate.evaluate_constraints_jacobian(problem);

   /******************************/
   /* build the symmetric matrix */
   /******************************/
   matrix.reset();

   // identity block
   for (size_t i = 0; i < number_variables; i++) {
      matrix.insert(1., i, i);
      matrix.finalize(i);
   }
   // Jacobian of general constraints
   for (size_t j = 0; j < problem.number_constraints; j++) {
      current_iterate.constraints_jacobian[j].for_each([&](size_t i, double derivative) {
         matrix.insert(derivative, i, number_variables + j);
      });
      matrix.finalize(number_variables + j);
   }
   DEBUG << "KKT matrix for least-square multipliers:\n" << matrix << "\n";

   /********************************/
   /* generate the right-hand side */
   /********************************/
   clear(rhs);

   // objective gradient
   current_iterate.objective_gradient.for_each([&](size_t i, double derivative) {
      rhs[i] += problem.objective_sign * derivative;
   });

   // variable bound constraints
   for (size_t i = 0; i < number_variables; i++) {
      rhs[i] -= current_iterate.multipliers.lower_bounds[i] + current_iterate.multipliers.upper_bounds[i];
   }

   // solve the system
   const size_t dimension = number_variables + problem.number_constraints;
   std::vector<double> solution(matrix.dimension);
   solver.factorize(dimension, matrix);
   solver.solve(dimension, matrix, rhs, solution);
   DEBUG << "Solution: "; print_vector(DEBUG, solution);

   // if least-square multipliers too big, discard them. Otherwise, store them
   if (norm_inf(solution, number_variables, problem.number_constraints) <= multipliers_max_size) {
      for (size_t j = 0; j < problem.number_constraints; j++) {
         multipliers[j] = solution[number_variables + j];
      }
   }
}