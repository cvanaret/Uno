#include <cmath>
#include <cassert>
#include "BQPDSolver.hpp"
#include "tools/Logger.hpp"
#include "linear_algebra/Vector.hpp"

#define BIG 1e30

extern "C" {
// fortran common block used in bqpd/bqpd.f
extern struct {
   int kk, ll, kkk, lll, mxws, mxlws;
} wsc_;

// fortran common for inertia correction in wdotd
extern struct {
   double alpha;
} kktalphac_;

extern void
bqpd_(const int* n, const int* m, int* k, int* kmax, double* a, int* la, double* x, double* bl, double* bu, double* f, double* fmin, double* g,
      double* r, double* w, double* e, int* ls, double* alp, int* lp, int* mlp, int* peq, double* ws, int* lws, const int* mode, int* ifail,
      int* info, int* iprint, int* nout);
}

// preallocate a bunch of stuff
BQPDSolver::BQPDSolver(size_t number_variables, size_t number_constraints, size_t maximum_number_nonzeros, bool quadratic_programming)
      : QPSolver(), number_variables(number_variables), number_constraints(number_constraints), maximum_number_nonzeros(maximum_number_nonzeros),
      lb(number_variables + number_constraints),
      ub(number_variables + number_constraints), jacobian(number_variables * (number_constraints + 1)),
      jacobian_sparsity(number_variables * (number_constraints + 1) + number_constraints + 3),
      kmax(quadratic_programming ? 500 : 0), alp(mlp), lp(mlp), ls(number_variables + number_constraints),
      w(number_variables + number_constraints), gradient_solution(number_variables), residuals(number_variables + number_constraints),
      e(number_variables + number_constraints),
      size_hessian_sparsity(quadratic_programming ? maximum_number_nonzeros + number_variables + 3 : 0),
      size_hessian_workspace(maximum_number_nonzeros + kmax * (kmax + 9) / 2 + 2 * number_variables + number_constraints + mxwk0),
      size_hessian_sparsity_workspace(size_hessian_sparsity + kmax + mxiwk0),
      hessian_values(size_hessian_workspace),
      hessian_sparsity(size_hessian_sparsity_workspace) {
   // active set
   for (size_t i = 0; i < this->number_variables + this->number_constraints; i++) {
      this->ls[i] = static_cast<int>(i + this->fortran_shift);
   }
}

Direction BQPDSolver::solve_QP(const std::vector<Range>& variables_bounds, const std::vector<Range>& constraints_bounds, const SparseVector<double>&
linear_objective, const std::vector<SparseVector<double>>& constraint_jacobian, const CSCSymmetricMatrix& hessian,
      const std::vector<double>& initial_point) {
   // Hessian
   for (size_t i = 0; i < hessian.number_nonzeros; i++) {
      this->hessian_values[i] = hessian.matrix[i];
   }
   // Hessian sparsity
   this->hessian_sparsity[0] = static_cast<int>(hessian.number_nonzeros + 1);
   for (size_t i = 0; i < hessian.number_nonzeros; i++) {
      this->hessian_sparsity[i + 1] = hessian.row_index[i] + static_cast<int>(this->fortran_shift);
   }
   for (size_t i = 0; i < hessian.dimension + 1; i++) {
      this->hessian_sparsity[hessian.number_nonzeros + i + 1] = hessian.column_start[i] + static_cast<int>(this->fortran_shift);
   }

   // if extra variables have been introduced, correct hessian.column_start
   // TODO move to HessianEvaluation
   size_t i = hessian.number_nonzeros + hessian.dimension + 2;
   const int last_value = hessian.column_start[hessian.dimension];
   for (size_t j = hessian.dimension; j < this->number_variables; j++) {
      this->hessian_sparsity[i] = last_value + static_cast<int>(this->fortran_shift);
      i++;
   }

   DEBUG << "Hessian: " << hessian;
   return this->solve_subproblem(variables_bounds, constraints_bounds, linear_objective, constraint_jacobian, initial_point);
}

Direction BQPDSolver::solve_LP(const std::vector<Range>& variables_bounds, const std::vector<Range>& constraints_bounds, const SparseVector<double>&
linear_objective, const std::vector<SparseVector<double>>& constraint_jacobian, const std::vector<double>& initial_point) {
   return this->solve_subproblem(variables_bounds, constraints_bounds, linear_objective, constraint_jacobian, initial_point);
}

bool all_variables_bounded(const std::vector<Range>& variables_bounds) {
   return std::all_of(variables_bounds.cbegin(), variables_bounds.cend(), [](const Range& bounds) {
      return -std::numeric_limits<double>::infinity() < bounds.lb && bounds.ub < std::numeric_limits<double>::infinity();
   });
}

void check_unboundedness(Status status, const std::vector<Range>& variables_bounds) {
   if (status == UNBOUNDED_PROBLEM) {
      if (all_variables_bounded(variables_bounds)) {
         assert(false && "BQPD: the problem is unbounded but the variables are bounded. A very negative objective value was probably reached");
      }
      else {
         assert(false && "BQPD: the problem is unbounded");
      }
   }
}

Direction BQPDSolver::solve_subproblem(const std::vector<Range>& variables_bounds, const std::vector<Range>& constraints_bounds,
      const SparseVector<double>& linear_objective, const std::vector<SparseVector<double>>& constraint_jacobian,
      const std::vector<double>& initial_point) {
   // initialize wsc_ common block (Hessian & workspace for bqpd)
   // setting the common block here ensures that several instances of BQPD can run simultaneously
   wsc_.kk = static_cast<int>(this->maximum_number_nonzeros);
   wsc_.ll = static_cast<int>(this->size_hessian_sparsity);
   wsc_.mxws = static_cast<int>(this->size_hessian_workspace);
   wsc_.mxlws = static_cast<int>(this->size_hessian_sparsity_workspace);
   kktalphac_.alpha = 0; // inertia control

   DEBUG << "objective gradient: ";
   DEBUG << linear_objective;
   for (size_t j = 0; j < constraint_jacobian.size(); j++) {
      DEBUG << "gradient c" << j << ": " << constraint_jacobian[j];
   }
   for (size_t i = 0; i < variables_bounds.size(); i++) {
      DEBUG << "Δx" << i << " in [" << variables_bounds[i].lb << ", " << variables_bounds[i].ub << "]\n";
   }
   for (size_t j = 0; j < constraints_bounds.size(); j++) {
      DEBUG << "linearized c" << j << " in [" << constraints_bounds[j].lb << ", " << constraints_bounds[j].ub << "]\n";
   }

   // Jacobian
   size_t current_index = 0;
   linear_objective.for_each([&](size_t i, double derivative) {
      this->jacobian[current_index] = derivative;
      this->jacobian_sparsity[current_index + 1] = static_cast<int>(i + this->fortran_shift);
      current_index++;
   });
   for (size_t j = 0; j < this->number_constraints; j++) {
      constraint_jacobian[j].for_each([&](size_t i, double derivative) {
         this->jacobian[current_index] = derivative;
         this->jacobian_sparsity[current_index + 1] = static_cast<int>(i + this->fortran_shift);
         current_index++;
      });
   }
   current_index++;
   this->jacobian_sparsity[0] = static_cast<int>(current_index);
   // header
   size_t size = 1;
   this->jacobian_sparsity[current_index] = static_cast<int>(size);
   current_index++;
   size += linear_objective.size();
   this->jacobian_sparsity[current_index] = static_cast<int>(size);
   current_index++;
   for (size_t j = 0; j < this->number_constraints; j++) {
      size += constraint_jacobian[j].size();
      this->jacobian_sparsity[current_index] = static_cast<int>(size);
      current_index++;
   }

   // bounds
   for (size_t i = 0; i < this->number_variables; i++) {
      this->lb[i] = (variables_bounds[i].lb == -std::numeric_limits<double>::infinity()) ? -BIG : variables_bounds[i].lb;
      this->ub[i] = (variables_bounds[i].ub == std::numeric_limits<double>::infinity()) ? BIG : variables_bounds[i].ub;
   }
   for (size_t j = 0; j < this->number_constraints; j++) {
      this->lb[this->number_variables + j] = (constraints_bounds[j].lb == -std::numeric_limits<double>::infinity()) ? -BIG : constraints_bounds[j].lb;
      this->ub[this->number_variables + j] = (constraints_bounds[j].ub == std::numeric_limits<double>::infinity()) ? BIG : constraints_bounds[j].ub;
   }

   Direction direction(this->number_variables, this->number_constraints);
   copy_from(direction.x, initial_point);
   const int n = static_cast<int>(this->number_variables);
   const int m = static_cast<int>(this->number_constraints);
   const int mode = static_cast<int>(this->mode);
   
   // solve the LP/QP
   bqpd_(&n, &m, &this->k, &this->kmax, this->jacobian.data(), this->jacobian_sparsity.data(),
         direction.x.data(), this->lb.data(), this->ub.data(), &direction.objective, &this->fmin, this->gradient_solution.data(),
         this->residuals.data(), this->w.data(), this->e.data(), this->ls.data(), this->alp.data(), this->lp.data(),
         &this->mlp, &this->peq_solution, this->hessian_values.data(), this->hessian_sparsity.data(), &mode, &this->ifail,
         this->info.data(), &this->iprint, &this->nout);
   direction.status = BQPDSolver::int_to_status(this->ifail);
   check_unboundedness(direction.status, variables_bounds);

   // project solution into bounds
   for (size_t i = 0; i < direction.x.size(); i++) {
      direction.x[i] = std::min(std::max(direction.x[i], variables_bounds[i].lb), variables_bounds[i].ub);
   }
   direction.norm = norm_inf(direction.x);
   this->analyze_constraints(direction);
   return direction;
}

void BQPDSolver::analyze_constraints(Direction& direction) {
   ConstraintPartition constraint_partition(this->number_constraints);

   // active constraints
   for (size_t j = 0; j < this->number_variables - this->k; j++) {
      size_t index = std::abs(this->ls[j]) - this->fortran_shift;

      if (index < this->number_variables) {
         // bound constraint
         if (0 <= this->ls[j]) { // lower bound active
            direction.multipliers.lower_bounds[index] = this->residuals[index];
            direction.active_set.bounds.at_lower_bound.push_back(index);
         }
         else { // upper bound active */
            direction.multipliers.upper_bounds[index] = -this->residuals[index];
            direction.active_set.bounds.at_upper_bound.push_back(index);
         }
      }
      else {
         // general constraint
         size_t constraint_index = index - this->number_variables;
         constraint_partition.feasible.push_back(constraint_index);
         constraint_partition.constraint_feasibility[constraint_index] = FEASIBLE;
         if (0 <= this->ls[j]) { // lower bound active
            direction.multipliers.constraints[constraint_index] = this->residuals[index];
            direction.active_set.constraints.at_lower_bound.push_back(constraint_index);
         }
         else { // upper bound active
            direction.multipliers.constraints[constraint_index] = -this->residuals[index];
            direction.active_set.constraints.at_upper_bound.push_back(constraint_index);
         }
      }
   }

   // inactive constraints
   for (size_t j = this->number_variables - this->k; j < this->number_variables + this->number_constraints; j++) {
      size_t index = std::abs(this->ls[j]) - this->fortran_shift;

      if (this->number_variables <= index) { // general constraints
         size_t constraint_index = index - this->number_variables;
         if (this->residuals[index] < 0.) { // infeasible constraint
            constraint_partition.infeasible.push_back(constraint_index);
            if (this->ls[j] < 0) { // upper bound violated
               constraint_partition.constraint_feasibility[constraint_index] = INFEASIBLE_UPPER;
            }
            else { // lower bound violated
               constraint_partition.constraint_feasibility[constraint_index] = INFEASIBLE_LOWER;
            }
         }
         else { // feasible constraint
            constraint_partition.feasible.push_back(constraint_index);
            constraint_partition.constraint_feasibility[constraint_index] = FEASIBLE;
         }
      }
   }
   direction.constraint_partition = constraint_partition;
}

Status BQPDSolver::int_to_status(int ifail) {
   assert(0 <= ifail && ifail <= 9 && "BQPDSolver.int_to_status: ifail does not belong to [0, 9]");
   return static_cast<Status> (ifail);
}