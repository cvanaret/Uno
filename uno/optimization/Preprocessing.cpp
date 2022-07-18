// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project root for details.

#include "Preprocessing.hpp"
#include "solvers/QP/BQPDSolver.hpp"
#include "tools/Range.hpp"

/*
void Preprocessing::enforce_linear_constraints(const Model& model, Iterate& first_iterate) {
   INFO << "Preprocessing phase: the problem has " << model.linear_constraints.size() << " linear constraints\n";
   if (!model.linear_constraints.empty()) {
      // count the infeasible constraints
      first_iterate.evaluate_constraints(model);
      int infeasible_linear_constraints = 0;
      model.linear_constraints.for_each_index([&](size_t j) {
         if (first_iterate.original_evaluations.constraints[j] < model.get_constraint_lower_bound(j) ||
             model.get_constraint_upper_bound(j) < first_iterate.original_evaluations.constraints[j]) {
            infeasible_linear_constraints++;
         }
      });
      INFO << "There are " << infeasible_linear_constraints << " infeasible linear constraints at the initial point\n";

      if (0 < infeasible_linear_constraints) {
         INFO << "Current point: "; print_vector(INFO, first_iterate.primals);
         const size_t number_linear_constraints = model.linear_constraints.size();
         BQPDSolver solver(model.number_variables, number_linear_constraints, model.number_variables, true);

         // objective: use a proximal term
         CSCSymmetricMatrix hessian = CSCSymmetricMatrix::identity(model.number_variables);
         SparseVector<double> linear_objective(0); // empty

         // constraints Jacobian
         std::vector<SparseVector<double>> constraint_jacobian(number_linear_constraints);
         for (size_t j = 0; j < model.number_constraints; j++) {
            constraint_jacobian[j].reserve(model.number_variables);
         }
         model.linear_constraints.for_each([&](size_t j, size_t linear_constraint_index) {
            model.evaluate_constraint_gradient(first_iterate.primals, j, constraint_jacobian[linear_constraint_index]);
         });

         // variables bounds
         std::vector<Interval> variables_bounds(model.number_variables);
         for (size_t i = 0; i < model.number_variables; i++) {
            variables_bounds[i] = {model.get_variable_lower_bound(i) - first_iterate.primals[i],
                  model.get_variable_upper_bound(i) - first_iterate.primals[i]};
         }

         // constraints bounds
         std::vector<Interval> constraint_bounds(number_linear_constraints);
         model.linear_constraints.for_each([&](size_t j, size_t linear_constraint_index) {
            constraint_bounds[linear_constraint_index] = {model.get_constraint_lower_bound(j) - first_iterate.original_evaluations.constraints[j],
                  model.get_constraint_upper_bound(j) - first_iterate.original_evaluations.constraints[j]};
         });

         // solve the convex QP
         std::vector<double> d0(model.number_variables);
         Direction direction = solver.solve_QP(model.number_variables, model.number_constraints, variables_bounds, constraint_bounds,
               linear_objective, constraint_jacobian, hessian, d0);
         if (direction.status == INFEASIBLE) {
            throw std::runtime_error("Linear constraints cannot be satisfied");
         }

         add_vectors(first_iterate.primals, direction.primals, 1., first_iterate.primals);
         // copy bound multipliers
         first_iterate.multipliers.lower_bounds = direction.multipliers.lower_bounds;
         first_iterate.multipliers.upper_bounds = direction.multipliers.upper_bounds;
         // copy constraint multipliers
         model.linear_constraints.for_each([&](size_t j, size_t linear_constraint_index) {
            first_iterate.multipliers.constraints[j] = direction.multipliers.constraints[linear_constraint_index];
         });
         INFO << "Linear feasible initial point: "; print_vector(INFO, first_iterate.primals); INFO << '\n';
      }
   }
}
 */

// compute a least-square approximation of the multipliers by solving a linear system (uses existing linear system)
void Preprocessing::compute_least_square_multipliers(const Model& model, SymmetricMatrix& matrix, std::vector<double>& rhs, LinearSolver& solver,
      Iterate& current_iterate, std::vector<double>& multipliers, double multipliers_max_norm) {
   const size_t number_variables = current_iterate.primals.size();
   current_iterate.evaluate_objective_gradient(model);
   current_iterate.evaluate_constraint_jacobian(model);

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
   for (size_t j = 0; j < model.number_constraints; j++) {
      current_iterate.original_evaluations.constraint_jacobian[j].for_each([&](size_t i, double derivative) {
         matrix.insert(derivative, i, number_variables + j);
      });
      matrix.finalize(number_variables + j);
   }
   DEBUG << "KKT matrix for least-square multipliers:\n" << matrix << '\n';

   /********************************/
   /* generate the right-hand side */
   /********************************/
   initialize_vector(rhs, 0.);
   // objective gradient
   current_iterate.original_evaluations.objective_gradient.for_each([&](size_t i, double derivative) {
      rhs[i] += model.objective_sign * derivative;
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
   DEBUG << "Solution: "; print_vector(DEBUG, solution, 0, matrix.dimension); DEBUG << '\n';

   // if least-square multipliers too big, discard them. Otherwise, store them
   if (norm_inf(solution, Range(number_variables, number_variables + model.number_constraints)) <= multipliers_max_norm) {
      for (size_t j = 0; j < model.number_constraints; j++) {
         multipliers[j] = solution[number_variables + j];
      }
   }
}