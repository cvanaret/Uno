// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "BQPDSolver.hpp"
#include "ingredients/local_models/LagrangeNewtonSubproblem.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Direction.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "reformulation/OptimizationProblem.hpp"
#include "tools/Infinity.hpp"
#include "tools/Logger.hpp"
#include "options/Options.hpp"
#include "fortran_interface.h"

#define WSC FC_GLOBAL(wsc,WSC)
#define ALPHAC FC_GLOBAL(alphac,ALPHAC)
#define BQPD FC_GLOBAL(bqpd,BQPD)
#define hessian_vector_product FC_GLOBAL(gdotx,GDOTX)

extern "C" {
   void hessian_vector_product(int *n, const double x[], const double ws[], const int lws[], double v[]);

   // fortran common block used in bqpd/bqpd.f
   extern struct {
      int kk, ll, kkk, lll, mxws, mxlws;
   } WSC;

   // fortran common for inertia correction in wdotd
   extern struct {
      double alpha;
   } ALPHAC;

   extern void
   BQPD(const int* n, const int* m, int* k, int* kmax, double* a, int* la, double* x, double* bl, double* bu, double* f, double* fmin, double* g,
         double* r, double* w, double* e, int* ls, double* alp, int* lp, int* mlp, int* peq, double* ws, int* lws, const int* mode, int* ifail,
         int* info, int* iprint, int* nout);
}

namespace uno {
   #define BIG 1e30

   // preallocate a bunch of stuff
   BQPDSolver::BQPDSolver(size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros, size_t number_jacobian_nonzeros,
         size_t number_hessian_nonzeros, BQPDProblemType problem_type, const Options& options):
         QPSolver(), number_hessian_nonzeros(number_hessian_nonzeros),
         lower_bounds(number_variables + number_constraints),
         upper_bounds(number_variables + number_constraints),
         jacobian(number_jacobian_nonzeros + number_objective_gradient_nonzeros), // Jacobian + objective gradient
         jacobian_sparsity(number_jacobian_nonzeros + number_objective_gradient_nonzeros + number_constraints + 3),
         kmax(problem_type == BQPDProblemType::QP ? options.get_int("BQPD_kmax") : 0), alp(static_cast<size_t>(this->mlp)),
         lp(static_cast<size_t>(this->mlp)),
         active_set(number_variables + number_constraints),
         w(number_variables + number_constraints), gradient_solution(number_variables), residuals(number_variables + number_constraints),
         e(number_variables + number_constraints),
         size_hessian_sparsity(problem_type == BQPDProblemType::QP ? number_hessian_nonzeros + number_variables + 3 : 0),
         size_hessian_workspace(number_hessian_nonzeros + static_cast<size_t>(this->kmax * (this->kmax + 9) / 2) + 2 * number_variables +
                                number_constraints + this->mxwk0),
         size_hessian_sparsity_workspace(this->size_hessian_sparsity + static_cast<size_t>(this->kmax) + this->mxiwk0),
         workspace(this->size_hessian_workspace),
         workspace_sparsity(this->size_hessian_sparsity_workspace),
         print_subproblem(options.get_bool("print_subproblem")),
         linear_objective(number_variables),
         constraints(number_constraints),
         constraint_jacobian(number_constraints, number_variables) {
      // default active set
      for (size_t variable_index: Range(number_variables + number_constraints)) {
         this->active_set[variable_index] = static_cast<int>(variable_index) + this->fortran_shift;
      }
   }

   void BQPDSolver::solve_QP(const LagrangeNewtonSubproblem& subproblem, const Vector<double>& initial_point, Direction& direction,
         const WarmstartInformation& warmstart_information) {
      if (this->print_subproblem) {
         DEBUG << "QP:\n";
      }
      this->set_up_subproblem(warmstart_information, subproblem);
      if (warmstart_information.objective_changed || warmstart_information.constraints_changed || warmstart_information.problem_changed) {
         this->save_hessian(subproblem);
      }
      this->solve_subproblem(initial_point, direction, warmstart_information, subproblem);
   }

   void BQPDSolver::solve_LP(const LagrangeNewtonSubproblem& subproblem, const Vector<double>& initial_point, Direction& direction,
         const WarmstartInformation& warmstart_information) {
      if (this->print_subproblem) {
         DEBUG << "LP:\n";
      }
      this->set_up_subproblem(warmstart_information, subproblem);
      this->solve_subproblem(initial_point, direction, warmstart_information, subproblem);
   }

   void BQPDSolver::set_up_subproblem(const WarmstartInformation& warmstart_information, const LagrangeNewtonSubproblem& subproblem) {
      // initialize WSC common block (Hessian & workspace for BQPD)
      // setting the common block here ensures that several instances of BQPD can run simultaneously
      WSC.kk = 0;
      WSC.ll = 0;
      WSC.mxws = static_cast<int>(this->size_hessian_workspace);
      WSC.mxlws = static_cast<int>(this->size_hessian_sparsity_workspace);
      ALPHAC.alpha = 0; // inertia control

      // objective gradient, constraints and constraint Jacobian
      if (warmstart_information.objective_changed) {
         subproblem.evaluate_objective_gradient(this->linear_objective);
      }
      if (warmstart_information.constraints_changed) {
         subproblem.evaluate_constraints(this->constraints);
         subproblem.evaluate_constraint_jacobian(this->constraint_jacobian);
      }
      if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {
         this->save_gradients_to_local_format(subproblem.number_constraints, this->linear_objective, this->constraint_jacobian);
      }

      // set bounds of the variable displacements
      if (warmstart_information.variable_bounds_changed) {
         subproblem.set_direction_bounds(this->lower_bounds, this->upper_bounds);
      }

      // set bounds of the linearized constraints
      if (warmstart_information.constraint_bounds_changed) {
         auto constraint_lower_bounds = view(this->lower_bounds, subproblem.number_variables, subproblem.number_variables + subproblem.number_constraints);
         auto constraint_upper_bounds = view(this->upper_bounds, subproblem.number_variables, subproblem.number_variables + subproblem.number_constraints);
         subproblem.set_constraint_bounds(this->constraints, constraint_lower_bounds, constraint_upper_bounds);
      }

      if (this->print_subproblem) {
         DEBUG << "objective gradient: " << linear_objective;
         for (size_t constraint_index: Range(subproblem.number_constraints)) {
            DEBUG << "gradient c" << constraint_index << ": " << constraint_jacobian[constraint_index];
         }

         for (size_t variable_index: Range(subproblem.number_variables)) {
            DEBUG << "d" << variable_index << " in [" << this->lower_bounds[variable_index] << ", " << this->upper_bounds[variable_index] << "]\n";
         }

         auto constraint_lower_bounds = view(this->lower_bounds, subproblem.number_variables, subproblem.number_variables + subproblem.number_constraints);
         auto constraint_upper_bounds = view(this->upper_bounds, subproblem.number_variables, subproblem.number_variables + subproblem.number_constraints);
         for (size_t constraint_index: Range(subproblem.number_constraints)) {
            DEBUG << "linearized c" << constraint_index << " in [" << constraint_lower_bounds[constraint_index] << ", " << constraint_upper_bounds[constraint_index] << "]\n";
         }
      }
   }

   void BQPDSolver::save_hessian(const LagrangeNewtonSubproblem& subproblem) {
      // hide subproblem in lws
      intptr_t pointer_to_subproblem = reinterpret_cast<intptr_t>(&subproblem);
      std::copy(reinterpret_cast<const char *>(&pointer_to_subproblem), reinterpret_cast<const char *>(&pointer_to_subproblem) + sizeof(intptr_t),
            reinterpret_cast<char *>(this->workspace_sparsity.data()));
      WSC.ll += sizeof(intptr_t);
   }

   void BQPDSolver::solve_subproblem(const Vector<double>& initial_point, Direction& direction,
         const WarmstartInformation& warmstart_information, const LagrangeNewtonSubproblem& subproblem) {
      direction.primals = initial_point;
      const int n = static_cast<int>(subproblem.number_variables);
      const int m = static_cast<int>(subproblem.number_constraints);

      const BQPDMode mode = this->determine_mode(warmstart_information);
      const int mode_integer = static_cast<int>(mode);

      // solve the LP/QP
      BQPD(&n, &m, &this->k, &this->kmax, this->jacobian.data(), this->jacobian_sparsity.data(), direction.primals.data(), this->lower_bounds.data(),
            this->upper_bounds.data(), &direction.subproblem_objective, &this->fmin, this->gradient_solution.data(), this->residuals.data(), this->w.data(),
            this->e.data(), this->active_set.data(), this->alp.data(), this->lp.data(), &this->mlp, &this->peq_solution, this->workspace.data(),
            this->workspace_sparsity.data(), &mode_integer, &this->ifail, this->info.data(), &this->iprint, &this->nout);
      const BQPDStatus bqpd_status = BQPDSolver::bqpd_status_from_int(this->ifail);
      direction.status = BQPDSolver::status_from_bqpd_status(bqpd_status);
      this->number_calls++;

      // project solution into bounds
      for (size_t variable_index: Range(subproblem.number_variables)) {
         direction.primals[variable_index] = std::min(std::max(direction.primals[variable_index], this->lower_bounds[variable_index]),
               this->upper_bounds[variable_index]);
      }
      this->categorize_constraints(subproblem.number_variables, direction);
   }

   BQPDMode BQPDSolver::determine_mode(const WarmstartInformation& warmstart_information) const {
      BQPDMode mode = (this->number_calls == 0) ? BQPDMode::ACTIVE_SET_EQUALITIES : BQPDMode::USER_DEFINED;
      // if problem changed, use cold start
      if (warmstart_information.problem_changed) {
         mode = BQPDMode::ACTIVE_SET_EQUALITIES;
      }
      // if only the variable bounds changed, reuse the active set estimate and the Jacobian information
      else if (warmstart_information.variable_bounds_changed && not warmstart_information.objective_changed &&
               not warmstart_information.constraints_changed && not warmstart_information.constraint_bounds_changed) {
         mode = BQPDMode::UNCHANGED_ACTIVE_SET_AND_JACOBIAN;
      }
      return mode;
   }

   void BQPDSolver::save_gradients_to_local_format(size_t number_constraints, const SparseVector<double>& linear_objective,
         const RectangularMatrix<double>& constraint_jacobian) {
      size_t current_index = 0;
      for (const auto [variable_index, derivative]: linear_objective) {
         this->jacobian[current_index] = derivative;
         this->jacobian_sparsity[current_index + 1] = static_cast<int>(variable_index) + this->fortran_shift;
         current_index++;
      }
      for (size_t constraint_index: Range(number_constraints)) {
         for (const auto [variable_index, derivative]: constraint_jacobian[constraint_index]) {
            this->jacobian[current_index] = derivative;
            this->jacobian_sparsity[current_index + 1] = static_cast<int>(variable_index) + this->fortran_shift;
            current_index++;
         }
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
      for (size_t constraint_index: Range(number_constraints)) {
         size += constraint_jacobian[constraint_index].size();
         this->jacobian_sparsity[current_index] = static_cast<int>(size);
         current_index++;
      }
   }

   void BQPDSolver::categorize_constraints(size_t number_variables, Direction& direction) {
      direction.multipliers.reset();

      // active constraints
      for (size_t active_constraint_index: Range(number_variables - static_cast<size_t>(this->k))) {
         const size_t index = static_cast<size_t>(std::abs(this->active_set[active_constraint_index]) - this->fortran_shift);

         if (index < number_variables) {
            // bound constraint
            if (0 <= this->active_set[active_constraint_index]) { // lower bound active
               direction.multipliers.lower_bounds[index] = this->residuals[index];
               direction.active_bounds.at_lower_bound.emplace_back(index);
            }
            else { // upper bound active */
               direction.multipliers.upper_bounds[index] = -this->residuals[index];
               direction.active_bounds.at_upper_bound.emplace_back(index);
            }
         }
         else {
            // general constraint
            size_t constraint_index = index - number_variables;
            if (0 <= this->active_set[active_constraint_index]) { // lower bound active
               direction.multipliers.constraints[constraint_index] = this->residuals[index];
            }
            else { // upper bound active
               direction.multipliers.constraints[constraint_index] = -this->residuals[index];
            }
         }
      }
   }

   BQPDStatus BQPDSolver::bqpd_status_from_int(int ifail) {
      assert(0 <= ifail && ifail <= 9 && "BQPDSolver.bqpd_status_from_int: ifail does not belong to [0, 9]");
      return static_cast<BQPDStatus>(ifail);
   }

   SubproblemStatus BQPDSolver::status_from_bqpd_status(BQPDStatus bqpd_status) {
      switch (bqpd_status) {
         case BQPDStatus::OPTIMAL:
            return SubproblemStatus::OPTIMAL;
         case BQPDStatus::UNBOUNDED_PROBLEM:
            return SubproblemStatus::UNBOUNDED_PROBLEM;
         case BQPDStatus::BOUND_INCONSISTENCY:
            DEBUG << "BQPD error: bound inconsistency\n";
            return SubproblemStatus::ERROR;
         case BQPDStatus::INFEASIBLE:
            return SubproblemStatus::INFEASIBLE;
            // errors
         case BQPDStatus::INCORRECT_PARAMETER:
            DEBUG << "BQPD error: incorrect parameter\n";
            return SubproblemStatus::ERROR;
         case BQPDStatus::LP_INSUFFICIENT_SPACE:
            DEBUG << "BQPD error: LP insufficient space\n";
            return SubproblemStatus::ERROR;
         case BQPDStatus::HESSIAN_INSUFFICIENT_SPACE:
            DEBUG << "BQPD kmax too small, continue anyway\n";
            return SubproblemStatus::ERROR;
         case BQPDStatus::SPARSE_INSUFFICIENT_SPACE:
            DEBUG << "BQPD error: sparse insufficient space\n";
            return SubproblemStatus::ERROR;
         case BQPDStatus::MAX_RESTARTS_REACHED:
            DEBUG << "BQPD max restarts reached\n";
            return SubproblemStatus::ERROR;
         case BQPDStatus::UNDEFINED:
            DEBUG << "BQPD error: undefined\n";
            return SubproblemStatus::ERROR;
      }
      throw std::invalid_argument("The BQPD ifail is not consistent with the Uno status values");
   }
} // namespace

void hessian_vector_product(int *n, const double x[], const double /*ws*/[], const int lws[], double v[]) {
   for (int i = 0; i < *n; i++) {
      v[i] = 0.;
   }

   /*
   int footer_start = lws[0];
   for (int i = 0; i < *n; i++) {
      for (int k = lws[footer_start + i]; k < lws[footer_start + i + 1]; k++) {
         int j = lws[k] - 1;
         v[i] += ws[k-1]*x[j];
         if (j != i) {
            // off-diagonal term
            v[j] += ws[k-1]*x[i];
         }
      }
   }
    */

   // retrieve subproblem
   intptr_t pointer_to_subproblem;
   std::copy(reinterpret_cast<const char *>(lws), reinterpret_cast<const char *>(lws) + sizeof(intptr_t), reinterpret_cast<char *>(&pointer_to_subproblem));
   const uno::LagrangeNewtonSubproblem* subproblem = reinterpret_cast<const uno::LagrangeNewtonSubproblem*>(pointer_to_subproblem);

   uno::Vector<double> my_x(*n);
   for (size_t variable_index: uno::Range(*n)) {
      my_x[variable_index] = x[variable_index];
   }
   uno::Vector<double> result(*n);
   subproblem->compute_hessian_vector_product(my_x, result);
   for (size_t variable_index: uno::Range(*n)) {
      v[variable_index] = result[variable_index];
   }
}