#include "Preprocessing.hpp"
#include "solvers/QP/BQPDSolver.hpp"

void Preprocessing::enforce_linear_constraints(const Problem& problem, Iterate& first_iterate) {
   INFO << "Preprocessing phase: the problem has " << problem.linear_constraints.size() << " linear constraints\n";
   if (!problem.linear_constraints.empty()) {
      // count the infeasible constraints
      first_iterate.evaluate_constraints(problem);
      int infeasible_linear_constraints = 0;
      problem.linear_constraints.for_each_index([&](size_t j) {
         if (first_iterate.constraints[j] < problem.get_constraint_lower_bound(j) ||
            problem.get_constraint_upper_bound(j) < first_iterate.constraints[j]) {
            infeasible_linear_constraints++;
         }
      });
      INFO << "There are " << infeasible_linear_constraints << " infeasible linear constraints at the initial point\n";

      if (0 < infeasible_linear_constraints) {
         INFO << "Current point: "; print_vector(INFO, first_iterate.x);
         const size_t number_constraints = problem.linear_constraints.size();
         BQPDSolver solver(problem.number_variables, number_constraints, problem.number_variables, true);

         // objective: use a proximal term
         CSCSymmetricMatrix hessian = CSCSymmetricMatrix::identity(problem.number_variables);
         SparseVector<double> linear_objective(0); // empty

         // constraints Jacobian
         std::vector<SparseVector<double>> constraint_jacobian(number_constraints);
         for (size_t j = 0; j < problem.number_constraints; j++) {
            constraint_jacobian[j].reserve(problem.number_variables);
         }
         problem.linear_constraints.for_each([&](size_t j, size_t linear_constraint_index) {
            problem.evaluate_constraint_gradient(first_iterate.x, j, constraint_jacobian[linear_constraint_index]);
         });

         // variables bounds
         std::vector<Range> variables_bounds(problem.number_variables);
         for (size_t i = 0; i < problem.number_variables; i++) {
            variables_bounds[i] = {problem.get_variable_lower_bound(i) - first_iterate.x[i],
                                   problem.get_variable_upper_bound(i) - first_iterate.x[i]};
         }

         // constraints bounds
         std::vector<Range> constraint_bounds(number_constraints);
         problem.linear_constraints.for_each([&](size_t j, size_t linear_constraint_index) {
            constraint_bounds[linear_constraint_index] = {problem.get_constraint_lower_bound(j) - first_iterate.constraints[j],
                                                           problem.get_constraint_upper_bound(j) - first_iterate.constraints[j]};
         });

         // solve the convex QP
         std::vector<double> d0(problem.number_variables);
         Direction direction = solver.solve_QP(variables_bounds, constraint_bounds, linear_objective, constraint_jacobian, hessian, d0);
         if (direction.status == INFEASIBLE) {
            throw std::runtime_error("Linear constraints cannot be satisfied");
         }

         add_vectors(first_iterate.x, direction.x, 1., first_iterate.x);
         // copy bound multipliers
         first_iterate.multipliers.lower_bounds = direction.multipliers.lower_bounds;
         first_iterate.multipliers.upper_bounds = direction.multipliers.upper_bounds;
         // copy constraint multipliers
         problem.linear_constraints.for_each([&](size_t j, size_t linear_constraint_index) {
            first_iterate.multipliers.constraints[j] = direction.multipliers.constraints[linear_constraint_index];
         });
         INFO << "Linear feasible initial point: "; print_vector(INFO, first_iterate.x); INFO << "\n";
      }
   }
}

// compute a least-square approximation of the multipliers by solving a linear system (uses existing linear system)
void Preprocessing::compute_least_square_multipliers(const Problem& problem, SymmetricMatrix& matrix,
      std::vector<double>& rhs, LinearSolver& solver, Iterate& current_iterate, std::vector<double>& multipliers, double multipliers_max_norm) {
   const size_t number_variables = current_iterate.x.size();
   const size_t number_constraints = multipliers.size();
   current_iterate.evaluate_objective_gradient(problem);
   current_iterate.evaluate_constraint_jacobian(problem);

   /******************************/
   /* build the symmetric matrix */
   /******************************/
   matrix.reset();
   matrix.dimension = number_variables + number_constraints;

   // identity block
   for (size_t i = 0; i < number_variables; i++) {
      matrix.insert(1., i, i);
      matrix.finalize(i);
   }
   // Jacobian of general constraints
   for (size_t j = 0; j < problem.number_constraints; j++) {
      current_iterate.constraint_jacobian[j].for_each([&](size_t i, double derivative) {
         matrix.insert(derivative, i, number_variables + j);
      });
      matrix.finalize(number_variables + j);
   }
   DEBUG << "KKT matrix for least-square multipliers:\n" << matrix << "\n";

   /********************************/
   /* generate the right-hand side */
   /********************************/
   initialize_vector(rhs, 0.);

   // objective gradient
   current_iterate.objective_gradient.for_each([&](size_t i, double derivative) {
      rhs[i] += problem.objective_sign * derivative;
   });

   // variable bound constraints
   for (size_t i = 0; i < number_variables; i++) {
      rhs[i] -= current_iterate.multipliers.lower_bounds[i] + current_iterate.multipliers.upper_bounds[i];
   }
   DEBUG << "RHS for least-square multipliers: "; print_vector(DEBUG, rhs, 0, matrix.dimension);

   /********************/
   /* solve the system */
   /********************/
   std::vector<double> solution(matrix.dimension);
   solver.factorize(matrix);
   solver.solve(matrix, rhs, solution);
   DEBUG << "Solution: "; print_vector(DEBUG, solution, 0, matrix.dimension);

   // if least-square multipliers too big, discard them. Otherwise, store them
   if (norm_inf(solution, number_variables, problem.number_constraints) <= multipliers_max_norm) {
      for (size_t j = 0; j < problem.number_constraints; j++) {
         multipliers[j] = solution[number_variables + j];
      }
   }
}