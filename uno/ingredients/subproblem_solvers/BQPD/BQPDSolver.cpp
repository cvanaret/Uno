// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "BQPDSolver.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/regularization_strategies/RegularizationStrategy.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/Logger.hpp"
#include "fortran_interface.h"

#define WSC FC_GLOBAL(wsc, WSC)
#define BQPD FC_GLOBAL(bqpd, BQPD)
#define hessian_vector_product FC_GLOBAL(gdotx, GDOTX)

extern "C" {
   void hessian_vector_product([[maybe_unused]] int *dimension, const double vector[], const double ws[],
      const int lws[], double result[]);

   // fortran common block used in bqpd/bqpd.f
   extern struct {
      int kk, ll, kkk, lll, mxws, mxlws;
   } WSC;

   extern void
   BQPD(const int* n, const int* m, int* k, int* kmax, double* a, int* la, double* x, double* bl, double* bu, double* f, double* fmin, double* g,
         double* r, double* w, double* e, int* ls, double* alp, int* lp, int* mlp, int* peq, double* ws, int* lws, const int* mode, int* ifail,
         int* info, int* iprint, int* nout);
}

namespace uno {
   #define BIG 1e30

   // heuristic to select kmax (the maximum size of the nullspace)
   // source: Minotaur code
   // https://github.com/coin-or/minotaur/blob/51a8bb78241a0b2e9ae94802aa76af8319f99192/src/interfaces/UnoEngine.cpp#L277
   int pick_kmax_heuristically(size_t number_variables, size_t number_constraints) {
      size_t kmax = 0;
      if (number_variables <= 1500) {
         kmax = number_variables;
      } else if (number_variables <= 5000) {
         kmax = number_variables - number_constraints/5;
      } else {
         kmax = 5000;
      }
      return static_cast<int>(kmax);
   }

   // preallocate a bunch of stuff
   BQPDSolver::BQPDSolver(const Options& options):
         QPSolver(),
         alp(static_cast<size_t>(this->mlp)),
         lp(static_cast<size_t>(this->mlp)),
         print_subproblem(options.get_bool("print_subproblem")) {
   }

   void BQPDSolver::initialize_memory(const OptimizationProblem& problem, const HessianModel& hessian_model,
         const RegularizationStrategy<double>& regularization_strategy) {
      this->w.resize(problem.number_variables + problem.number_constraints);
      this->gradient_solution.resize(problem.number_variables);
      this->residuals.resize(problem.number_variables + problem.number_constraints);
      this->e.resize(problem.number_variables + problem.number_constraints);

      this->lower_bounds.resize(problem.number_variables + problem.number_constraints);
      this->upper_bounds.resize(problem.number_variables + problem.number_constraints);
      this->constraints.resize(problem.number_constraints);
      this->linear_objective.reserve(problem.number_objective_gradient_nonzeros());
      this->constraint_jacobian.resize(problem.number_constraints, problem.number_variables);
      // Jacobian + objective gradient
      this->bqpd_jacobian.resize(problem.number_jacobian_nonzeros() + problem.number_objective_gradient_nonzeros());
      this->bqpd_jacobian_sparsity.resize(problem.number_jacobian_nonzeros() + problem.number_objective_gradient_nonzeros() +
         problem.number_constraints + 3);
      // default active set
      this->active_set.resize(problem.number_variables + problem.number_constraints);
      for (size_t variable_index: Range(problem.number_variables + problem.number_constraints)) {
         this->active_set[variable_index] = static_cast<int>(variable_index) + this->fortran_shift;
      }

      // Hessian: determine whether the subproblem has curvature
      const size_t number_hessian_nonzeros = problem.number_hessian_nonzeros(hessian_model);
      const size_t regularization_size = (!hessian_model.is_positive_definite() &&
         regularization_strategy.performs_primal_regularization()) ? problem.get_number_original_variables() : 0;
      // if the Hessian model only has an explicit representation, allocate an explicit Hessian matrix
      if (!hessian_model.has_implicit_representation() && hessian_model.has_explicit_representation()) {
         this->hessian = SparseSymmetricMatrix<COOFormat<size_t, double>>(problem.number_variables,
            number_hessian_nonzeros, regularization_size);
      }
      const size_t number_regularized_hessian_nonzeros = number_hessian_nonzeros + regularization_size;
      this->kmax = (0 < number_regularized_hessian_nonzeros) ? pick_kmax_heuristically(problem.number_variables,
         problem.number_constraints) : 0;

      // allocation of integer and real workspaces
      this->mxws = static_cast<size_t>(this->kmax * (this->kmax + 9) / 2) + 2 * problem.number_variables +
         problem.number_constraints /* (required by bqpd.f) */ + 5 * problem.number_variables + this->nprof /* (required
         by sparseL.f) */;
      // 4 pointers hidden in lws
      constexpr size_t hidden_pointers_size = 4*sizeof(intptr_t);
      this->mxlws = hidden_pointers_size + static_cast<size_t>(this->kmax) /* (required by bqpd.f) */ +
         9 * problem.number_variables + problem.number_constraints /* (required by sparseL.f) */;
      this->ws.resize(this->mxws);
      this->lws.resize(this->mxlws);
   }

   void BQPDSolver::solve_inequality_constrained_subproblem(Statistics& statistics, Subproblem& subproblem, const Vector<double>& initial_point,
         Direction& direction, const WarmstartInformation& warmstart_information) {
      this->set_up_subproblem(statistics, subproblem, warmstart_information);
      if (this->print_subproblem) {
         this->display_subproblem(subproblem, initial_point);
      }
      this->solve_subproblem(subproblem, initial_point, direction, warmstart_information);
   }

   double BQPDSolver::hessian_quadratic_product(const Vector<double>& vector) const {
      return this->hessian.quadratic_product(vector, vector);
   }

   // protected member functions

   void BQPDSolver::set_up_subproblem(Statistics& statistics, const Subproblem& subproblem,
         const WarmstartInformation& warmstart_information) {
      // initialize wsc_ common block (Hessian & workspace for BQPD)
      // setting the common block here ensures that several instances of BQPD can run simultaneously
      WSC.mxws = static_cast<int>(this->mxws);
      WSC.mxlws = static_cast<int>(this->mxlws);

      // evaluate the functions based on warmstart information
      if (warmstart_information.objective_changed) {
         subproblem.evaluate_objective_gradient(this->linear_objective);
      }
      if (warmstart_information.constraints_changed) {
         subproblem.evaluate_constraints(this->constraints);
         subproblem.evaluate_jacobian(this->constraint_jacobian);
      }
      if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {
         this->evaluate_hessian = true;
      }

      // variable bounds
      if (warmstart_information.variable_bounds_changed) {
         subproblem.set_variables_bounds(this->lower_bounds, this->upper_bounds);
      }

      // constraint bounds
      if (warmstart_information.constraint_bounds_changed || warmstart_information.constraints_changed) {
         auto constraints_lower_bounds = view(this->lower_bounds, subproblem.number_variables, subproblem.number_variables + subproblem.number_constraints);
         auto constraints_upper_bounds = view(this->upper_bounds, subproblem.number_variables, subproblem.number_variables + subproblem.number_constraints);
         subproblem.set_constraints_bounds(constraints_lower_bounds, constraints_upper_bounds, this->constraints);
      }

      // replace INFs with large finite values (TODO: is that really useful?)
      for (size_t variable_index: Range(subproblem.number_variables + subproblem.number_constraints)) {
         this->lower_bounds[variable_index] = std::max(-BIG, this->lower_bounds[variable_index]);
         this->upper_bounds[variable_index] = std::min(BIG, this->upper_bounds[variable_index]);
      }

      // save objective gradient and constraint Jacobian into BQPD workspace
      if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {
         this->save_gradients_into_workspace(subproblem.number_constraints);
      }
      this->hide_pointers_in_workspace(statistics, subproblem);
   }

   void BQPDSolver::display_subproblem(const Subproblem& subproblem, const Vector<double>& initial_point) const {
      DEBUG << "Subproblem:\n";
      DEBUG << "Hessian: " << this->hessian;
      DEBUG << "objective gradient: " << this->linear_objective;
      for (size_t constraint_index: Range(subproblem.number_constraints)) {
         DEBUG << "gradient c" << constraint_index << ": " << this->constraint_jacobian[constraint_index];
      }
      for (size_t variable_index: Range(subproblem.number_variables)) {
         DEBUG << "d" << variable_index << " in [" << this->lower_bounds[variable_index] << ", " << this->upper_bounds[variable_index] << "]\n";
      }
      for (size_t constraint_index: Range(subproblem.number_constraints)) {
         DEBUG << "linearized c" << constraint_index << " in [" << this->lower_bounds[subproblem.number_variables + constraint_index] << ", " <<
            this->upper_bounds[subproblem.number_variables + constraint_index] << "]\n";
      }
      DEBUG << "Initial point: " << initial_point << '\n';
   }

   void BQPDSolver::solve_subproblem(const Subproblem& subproblem, const Vector<double>& initial_point, Direction& direction,
         const WarmstartInformation& warmstart_information) {
      direction.primals = initial_point;
      const int n = static_cast<int>(subproblem.number_variables);
      const int m = static_cast<int>(subproblem.number_constraints);

      const BQPDMode mode = BQPDSolver::determine_mode(warmstart_information);
      const int mode_integer = static_cast<int>(mode);

      // solve the LP/QP
      DEBUG2 << "Running BQPD\n";
      BQPD(&n, &m, &this->k, &this->kmax, this->bqpd_jacobian.data(), this->bqpd_jacobian_sparsity.data(), direction.primals.data(),
         this->lower_bounds.data(), this->upper_bounds.data(), &direction.subproblem_objective, &this->fmin, this->gradient_solution.data(),
         this->residuals.data(), this->w.data(), this->e.data(), this->active_set.data(), this->alp.data(), this->lp.data(), &this->mlp,
         &this->peq_solution, this->ws.data(), this->lws.data(), &mode_integer, &this->ifail, this->info.data(),
         &this->iprint, &this->nout);
      DEBUG2 << "Ran BQPD\n";
      const BQPDStatus bqpd_status = BQPDSolver::bqpd_status_from_int(this->ifail);
      direction.status = BQPDSolver::status_from_bqpd_status(bqpd_status);

      // project solution into bounds
      for (size_t variable_index: Range(subproblem.number_variables)) {
         direction.primals[variable_index] = std::min(std::max(direction.primals[variable_index], this->lower_bounds[variable_index]),
            this->upper_bounds[variable_index]);
      }
      this->set_multipliers(subproblem.number_variables, direction.multipliers);
   }

   BQPDMode BQPDSolver::determine_mode(const WarmstartInformation& warmstart_information) {
      BQPDMode mode = BQPDMode::USER_DEFINED;
      // if problem structure changed, use cold start
      if (warmstart_information.hessian_sparsity_changed || warmstart_information.jacobian_sparsity_changed) {
         mode = BQPDMode::ACTIVE_SET_EQUALITIES;
      }
      // if only the variable bounds changed, reuse the active set estimate and the Jacobian information
      else if (warmstart_information.variable_bounds_changed && !warmstart_information.objective_changed &&
               !warmstart_information.constraints_changed && !warmstart_information.constraint_bounds_changed) {
         mode = BQPDMode::UNCHANGED_ACTIVE_SET_AND_JACOBIAN;
      }
      return mode;
   }

   // hide pointers to arbitrary objects into this->workspace_sparsity (BQPD's lws)
   void BQPDSolver::hide_pointers_in_workspace(Statistics& statistics, const Subproblem& subproblem) {
      WSC.kk = 0; // length of ws that is used by gdotx
      WSC.ll = 0; // length of lws that is used by gdotx

      // hide pointer to this->evaluate_hessian, statistics, subproblem and Hessian
      hide_pointer(0, this->lws.data(), this->evaluate_hessian);
      WSC.ll += sizeof(intptr_t);
      hide_pointer(1, this->lws.data(), statistics);
      WSC.ll += sizeof(intptr_t);
      hide_pointer(2, this->lws.data(), subproblem);
      WSC.ll += sizeof(intptr_t);
      hide_pointer(3, this->lws.data(), this->hessian);
      WSC.ll += sizeof(intptr_t);
   }

   void BQPDSolver::save_gradients_into_workspace(size_t number_constraints) {
      size_t current_index = 0;
      for (const auto [variable_index, derivative]: this->linear_objective) {
         assert(current_index < this->bqpd_jacobian.size() && "The allocation of bqpd_jacobian was not sufficient");
         assert(current_index + 1 < this->bqpd_jacobian_sparsity.size() && "The allocation of bqpd_jacobian_sparsity was not sufficient");
         this->bqpd_jacobian[current_index] = derivative;
         this->bqpd_jacobian_sparsity[current_index + 1] = static_cast<int>(variable_index) + this->fortran_shift;
         current_index++;
      }
      for (size_t constraint_index: Range(number_constraints)) {
         for (const auto [variable_index, derivative]: this->constraint_jacobian[constraint_index]) {
            assert(current_index < this->bqpd_jacobian.size() && "The allocation of bqpd_jacobian was not sufficient");
            assert(current_index + 1 < this->bqpd_jacobian_sparsity.size() && "The allocation of bqpd_jacobian_sparsity was not sufficient");
            this->bqpd_jacobian[current_index] = derivative;
            this->bqpd_jacobian_sparsity[current_index + 1] = static_cast<int>(variable_index) + this->fortran_shift;
            current_index++;
         }
      }
      current_index++;
      this->bqpd_jacobian_sparsity[0] = static_cast<int>(current_index);
      // header
      size_t size = 1;
      this->bqpd_jacobian_sparsity[current_index] = static_cast<int>(size);
      current_index++;
      size += this->linear_objective.size();
      this->bqpd_jacobian_sparsity[current_index] = static_cast<int>(size);
      current_index++;
      for (size_t constraint_index: Range(number_constraints)) {
         size += this->constraint_jacobian[constraint_index].size();
         this->bqpd_jacobian_sparsity[current_index] = static_cast<int>(size);
         current_index++;
      }
   }

   void BQPDSolver::set_multipliers(size_t number_variables, Multipliers& direction_multipliers) const {
      direction_multipliers.reset();
      // active constraints
      for (size_t active_constraint_index: Range(number_variables - static_cast<size_t>(this->k))) {
         const size_t index = static_cast<size_t>(std::abs(this->active_set[active_constraint_index]) - this->fortran_shift);

         if (index < number_variables) {
            // bound constraint
            if (0 <= this->active_set[active_constraint_index]) { // lower bound active
               direction_multipliers.lower_bounds[index] = this->residuals[index];
            }
            else { // upper bound active */
               direction_multipliers.upper_bounds[index] = -this->residuals[index];
            }
         }
         else {
            // general constraint
            size_t constraint_index = index - number_variables;
            if (0 <= this->active_set[active_constraint_index]) { // lower bound active
               direction_multipliers.constraints[constraint_index] = this->residuals[index];
            }
            else { // upper bound active
               direction_multipliers.constraints[constraint_index] = -this->residuals[index];
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

void hessian_vector_product(int* dimension, const double vector[], const double /*ws*/[], const int lws[], double result[]) {
   assert(dimension != nullptr && "BQPDSolver::hessian_vector_product: the dimension n passed by pointer is NULL");

   for (size_t i = 0; i < static_cast<size_t>(*dimension); i++) {
      result[i] = 0.;
   }

   // retrieve flag evaluate_hessian, statistics, subproblem and Hessian
   bool* evaluate_hessian = uno::retrieve_pointer<bool>(0, lws);
   uno::Statistics* statistics = uno::retrieve_pointer<uno::Statistics>(1, lws);
   uno::Subproblem* subproblem = uno::retrieve_pointer<uno::Subproblem>(2, lws);
   uno::SymmetricMatrix<size_t, double>* hessian = uno::retrieve_pointer<uno::SymmetricMatrix<size_t, double>>(3, lws);
   assert(evaluate_hessian != nullptr && "BQPD's hessian_vector_product: the flag evaluate_hessian is NULL");
   assert(statistics != nullptr && "BQPD's hessian_vector_product: statistics is NULL");
   assert(subproblem != nullptr && "BQPD's hessian_vector_product: subproblem is NULL");
   assert(hessian != nullptr && "BQPD's hessian_vector_product: hessian is NULL");

   // by default, try to perform a Hessian-vector product if possible
   if (subproblem->has_implicit_hessian_representation()) {
      subproblem->compute_hessian_vector_product(vector, result);
   }
   // otherwise, try to compute the explicit matrix
   else if (subproblem->has_explicit_hessian_representation()) {
      // if the Hessian has not been evaluated at the current point, evaluate it
      if (*evaluate_hessian) {
         hessian->reset();
         subproblem->compute_regularized_hessian(*statistics, *hessian);
         *evaluate_hessian = false;
      }
      hessian->product(vector, result);
   }
   else {
      throw std::runtime_error("Hessian model cannot be evaluated");
   }
}