#include <cmath>
#include <cassert>
#include "BQPDSolver.hpp"
#include "tools/Logger.hpp"
#include "linear_algebra/Vector.hpp"

#define BIG 1e30

extern "C" {
   /* fortran common block used in bqpd/bqpd.f */
   extern struct {
      int kk, ll, kkk, lll, mxws, mxlws;
   } wsc_;

   /* fortran common for inertia correction in wdotd */
   extern struct {
      double alpha;
   } kktalphac_;

   extern void
   bqpd_(const int* n, const int* m, int* k, int* kmax, double* a, int* la, double* x, double* bl, double* bu, double* f, double* fmin, double* g,
         double* r, double* w, double* e, int* ls, double* alp, int* lp, int* mlp, int* peq, double* ws, int* lws, const int* mode, int* ifail,
         int* info, int* iprint, int* nout);
}

/* preallocate a bunch of stuff */
BQPDSolver::BQPDSolver(size_t number_variables, size_t number_constraints, size_t maximum_number_nonzeros, bool quadratic_programming)
      : QPSolver(), n(number_variables), m(number_constraints), maximum_number_nonzeros(maximum_number_nonzeros), lb(n + m),
      ub(n + m), jacobian(n * (m + 1)), jacobian_sparsity(n * (m + 1) + m + 3),
      kmax(quadratic_programming ? 500 : 0), alp(mlp), lp(mlp), ls(n + m),
      w(n + m), gradient_solution(n), residuals(n + m), e(n + m),
      size_hessian_sparsity(quadratic_programming ? maximum_number_nonzeros + n + 3 : 0),
      size_hessian_workspace(maximum_number_nonzeros + kmax * (kmax + 9) / 2 + 2 * n + m + mxwk0),
      size_hessian_sparsity_workspace(size_hessian_sparsity + kmax + mxiwk0),
      hessian(size_hessian_workspace),
      hessian_sparsity(size_hessian_sparsity_workspace),
      solution(number_variables) {
   // active set
   for (size_t i = 0; i < this->n + this->m; i++) {
      this->ls[i] = static_cast<int>(i + this->fortran_shift);
   }
}

Direction BQPDSolver::solve_QP(const std::vector<Range>& variables_bounds, const std::vector<Range>& constraints_bounds, const SparseVector<double>&
   linear_objective, const std::vector<SparseVector<double>>& constraints_jacobian, const CSCSymmetricMatrix& hessian, const std::vector<double>& initial_point) {
   /* Hessian */
   for (size_t i = 0; i < hessian.number_nonzeros; i++) {
      this->hessian[i] = hessian.matrix[i];
   }
   /* Hessian sparsity */
   this->hessian_sparsity[0] = static_cast<int>(hessian.number_nonzeros + 1);
   for (size_t i = 0; i < hessian.number_nonzeros; i++) {
      this->hessian_sparsity[i + 1] = static_cast<int>(hessian.row_index[i] + this->fortran_shift);
   }
   for (size_t i = 0; i < hessian.dimension + 1; i++) {
      this->hessian_sparsity[hessian.number_nonzeros + i + 1] = static_cast<int>(hessian.column_start[i] + this->fortran_shift);
   }

   // if extra variables have been introduced, correct hessian.column_start
   size_t i = hessian.number_nonzeros + hessian.dimension + 2;
   size_t last_value = hessian.column_start[hessian.dimension];
   for (size_t j = hessian.dimension; j < this->n; j++) {
      this->hessian_sparsity[i] = static_cast<int>(last_value + this->fortran_shift);
      i++;
   }

   DEBUG << "Hessian: " << hessian;
   return this->solve_subproblem(variables_bounds, constraints_bounds, linear_objective, constraints_jacobian, initial_point);
}

Direction BQPDSolver::solve_LP(const std::vector<Range>& variables_bounds, const std::vector<Range>& constraints_bounds, const SparseVector<double>&
   linear_objective, const std::vector<SparseVector<double>>& constraints_jacobian, const std::vector<double>& initial_point) {
   return this->solve_subproblem(variables_bounds, constraints_bounds, linear_objective, constraints_jacobian, initial_point);
}

Direction BQPDSolver::solve_subproblem(const std::vector<Range>& variables_bounds, const std::vector<Range>& constraints_bounds,
      const SparseVector<double>& linear_objective, const std::vector<SparseVector<double>>& constraints_jacobian, const std::vector<double>& x) {
   /* initialize wsc_ common block (Hessian & workspace for bqpd) */
   // setting the common block here ensures that several instances of BQPD can run simultaneously
   wsc_.kk = static_cast<int>(this->maximum_number_nonzeros);
   wsc_.ll = static_cast<int>(this->size_hessian_sparsity);
   wsc_.mxws = static_cast<int>(this->size_hessian_workspace);
   wsc_.mxlws = static_cast<int>(this->size_hessian_sparsity_workspace);
   kktalphac_.alpha = 0; // inertia control

   DEBUG << "objective gradient: ";
   DEBUG << linear_objective;
   for (size_t j = 0; j < constraints_jacobian.size(); j++) {
      DEBUG << "gradient c" << j << ": " << constraints_jacobian[j];
   }
   for (size_t i = 0; i < variables_bounds.size(); i++) {
      DEBUG << "Δx" << i << " in [" << variables_bounds[i].lb << ", " << variables_bounds[i].ub << "]\n";
   }
   for (size_t j = 0; j < constraints_bounds.size(); j++) {
      DEBUG << "linearized c" << j << " in [" << constraints_bounds[j].lb << ", " << constraints_bounds[j].ub << "]\n";
   }

   /* Jacobian */
   int current_index = 0;
   linear_objective.for_each([&](size_t i, double derivative) {
      this->jacobian[current_index] = derivative;
      this->jacobian_sparsity[current_index + 1] = static_cast<int>(i + this->fortran_shift);
      current_index++;
   });
   for (size_t j = 0; j < this->m; j++) {
      constraints_jacobian[j].for_each([&](size_t i, double derivative) {
         this->jacobian[current_index] = derivative;
         this->jacobian_sparsity[current_index + 1] = static_cast<int>(i + this->fortran_shift);
         current_index++;
      });
   }
   current_index++;
   this->jacobian_sparsity[0] = current_index;
   // header
   size_t size = 1;
   this->jacobian_sparsity[current_index] = static_cast<int>(size);
   current_index++;
   size += linear_objective.size();
   this->jacobian_sparsity[current_index] = static_cast<int>(size);
   current_index++;
   for (size_t j = 0; j < this->m; j++) {
      size += constraints_jacobian[j].size();
      this->jacobian_sparsity[current_index] = static_cast<int>(size);
      current_index++;
   }

   /* bounds */
   for (size_t i = 0; i < this->n; i++) {
      this->lb[i] = (variables_bounds[i].lb == -INFINITY) ? -BIG : variables_bounds[i].lb;
      this->ub[i] = (variables_bounds[i].ub == INFINITY) ? BIG : variables_bounds[i].ub;
   }
   for (size_t j = 0; j < this->m; j++) {
      this->lb[this->n + j] = (constraints_bounds[j].lb == -INFINITY) ? -BIG : constraints_bounds[j].lb;
      this->ub[this->n + j] = (constraints_bounds[j].ub == INFINITY) ? BIG : constraints_bounds[j].ub;
   }

   /* initial point */
   copy_from(this->solution, x);

   /* call BQPD */
   const int n = static_cast<int>(this->n);
   const int m = static_cast<int>(this->m);
   const int mode = static_cast<int>(this->mode);
   bqpd_(&n, &m, &this->k, &this->kmax, this->jacobian.data(), this->jacobian_sparsity.data(),
         this->solution.data(),this->lb.data(), this->ub.data(), &this->f_solution, &this->fmin, this->gradient_solution.data(),
         this->residuals.data(),this->w.data(), this->e.data(), this->ls.data(), this->alp.data(), this->lp.data(),
         &this->mlp, &this->peq_solution,this->hessian.data(), this->hessian_sparsity.data(), &mode, &this->ifail,
         this->info.data(), &this->iprint, &this->nout);

   /* project solution into bounds */
   for (unsigned int i = 0; i < this->solution.size(); i++) {
      if (this->solution[i] < variables_bounds[i].lb) {
         this->solution[i] = variables_bounds[i].lb;
      }
      else if (variables_bounds[i].ub < this->solution[i]) {
         this->solution[i] = variables_bounds[i].ub;
      }
   }
   return this->generate_direction();
}

Direction BQPDSolver::generate_direction() {
   Direction direction(this->n, this->m);
   copy_from(direction.x, this->solution);
   direction.status = this->int_to_status(this->ifail);
   assert(direction.status != UNBOUNDED_PROBLEM && "BQPD: the problem is unbounded");
   direction.norm = norm_inf(this->solution);
   direction.objective = this->f_solution;

   /* active constraints */
   for (size_t j = 0; j < this->n - this->k; j++) {
      size_t index = std::abs(this->ls[j]) - this->fortran_shift;

      if (index < this->n) {
         // bound constraint
         if (0 <= this->ls[j]) { /* lower bound active */
            direction.multipliers.lower_bounds[index] = this->residuals[index];
            direction.active_set.bounds.at_lower_bound.push_back(index);
         }
         else { /* upper bound active */
            direction.multipliers.upper_bounds[index] = -this->residuals[index];
            direction.active_set.bounds.at_upper_bound.push_back(index);
         }
      }
      else {
         // general constraint
         size_t constraint_index = index - this->n;
         direction.constraint_partition.feasible.push_back(constraint_index);
         direction.constraint_partition.constraint_feasibility[constraint_index] = FEASIBLE;
         if (0 <= this->ls[j]) { /* lower bound active */
            direction.multipliers.constraints[constraint_index] = this->residuals[index];
            direction.active_set.constraints.at_lower_bound.push_back(constraint_index);
         }
         else { /* upper bound active */
            direction.multipliers.constraints[constraint_index] = -this->residuals[index];
            direction.active_set.constraints.at_upper_bound.push_back(constraint_index);
         }
      }
   }

   /* inactive constraints */
   for (size_t j = this->n - this->k; j < this->n + this->m; j++) {
      size_t index = std::abs(this->ls[j]) - this->fortran_shift;

      if (this->n <= index) { // general constraints
         size_t constraint_index = index - this->n;
         if (this->residuals[index] < 0.) { // infeasible constraint
            direction.constraint_partition.infeasible.push_back(constraint_index);
            if (this->ls[j] < 0) { // upper bound violated
               direction.constraint_partition.constraint_feasibility[constraint_index] = INFEASIBLE_UPPER;
            }
            else { // lower bound violated
               direction.constraint_partition.constraint_feasibility[constraint_index] = INFEASIBLE_LOWER;
            }
         }
         else { // feasible constraint
            direction.constraint_partition.feasible.push_back(constraint_index);
            direction.constraint_partition.constraint_feasibility[constraint_index] = FEASIBLE;
         }
      }
   }
   return direction;
}

Status BQPDSolver::int_to_status(int ifail) {
   assert(0 <= ifail && ifail <= 9 && "BQPDSolver.int_to_status: ifail does not belong to [0, 9]");
   return static_cast<Status> (ifail);
}