// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "BQPDSolver.hpp"
#include "BQPDQuadraticProgram.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/Indexing.hpp"
#include "linear_algebra/Vector.hpp"
#include "linear_algebra/VectorView.hpp"
#include "optimization/Direction.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "symbolic/Range.hpp"
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

   extern void BQPD(const int* n, const int* m, int* k, int* kmax, double* a, int* la, double* x, double* bl, double* bu,
      double* f, double* fmin, double* g, double* r, double* w, double* e, int* ls, double* alp, int* lp, int* mlp, int* peq,
      double* ws, int* lws, const int* mode, int* ifail, int* info, int* iprint, int* nout);
}

namespace uno {
   // heuristics to select kmax (the maximum size of the nullspace)

   // approach implemented in filterSQP
   int pick_kmax_filtersqp(size_t number_variables, size_t /*number_constraints*/) {
      return std::min(static_cast<int>(number_variables), 500);
   }

   // approach implemented in Minotaur
   // https://github.com/coin-or/minotaur/blob/51a8bb78241a0b2e9ae94802aa76af8319f99192/src/interfaces/UnoEngine.cpp#L277
   int pick_kmax_minotaur(size_t number_variables, size_t number_constraints) {
      size_t kmax = 0;
      if (number_variables <= 1500) {
         kmax = number_variables;
      }
      else if (number_variables <= 5000) {
         kmax = number_variables - number_constraints/5;
      }
      else {
         kmax = 5000;
      }
      return static_cast<int>(kmax);
   }

   // preallocate a bunch of stuff
   BQPDSolver::BQPDSolver(const Options& options):
         QPSolver(),
         // select a heuristic to pick kmax (the max size of the nullspace)
         pick_kmax(options.get_string("BQPD_kmax_heuristic") == "minotaur" ? pick_kmax_minotaur : pick_kmax_filtersqp),
         alp(static_cast<size_t>(this->mlp)),
         lp(static_cast<size_t>(this->mlp)),
         print_subproblem(options.get_bool("print_subproblem")) {
      // construct an empty BQPD-native quadratic program so that get_quadratic_program() can be used to build
      // it directly from data (no Subproblem); the full solver instead calls initialize_memory(subproblem)
      this->quadratic_program = std::make_unique<BQPDQuadraticProgram>();
   }

   BQPDSolver::~BQPDSolver() = default;

   void BQPDSolver::initialize_memory(const Subproblem& subproblem) {
      this->number_variables = subproblem.number_variables;
      this->number_constraints = subproblem.number_constraints;
      // build the BQPD-native quadratic program (allocates gradients/Hessian/bounds storage and the
      // iteration-invariant sparsity patterns)
      this->quadratic_program = std::make_unique<BQPDQuadraticProgram>();
      this->quadratic_program->initialize_memory(subproblem);
      // allocate the BQPD algorithm scratch
      this->allocate_workspace(subproblem.number_variables, subproblem.number_constraints, subproblem.number_jacobian_nonzeros(),
         subproblem.has_curvature());
   }

   void BQPDSolver::allocate_workspace(size_t number_variables, size_t number_constraints, size_t number_jacobian_nonzeros, bool is_qp) {
      // BQPD algorithm scratch
      this->w.resize(number_variables + number_constraints);
      this->gradient_solution.resize(number_variables);
      this->residuals.resize(number_variables + number_constraints);
      this->e.resize(number_variables + number_constraints);

      // default active set
      this->active_set.resize(number_variables + number_constraints);
      for (size_t variable_index: Range(number_variables + number_constraints)) {
         this->active_set[variable_index] = static_cast<int>(variable_index + Indexing::Fortran_indexing);
      }

      // kmax (the max size of the nullspace) is 0 for an LP
      this->kmax = is_qp ? pick_kmax(number_variables, number_constraints) : 0;

      // allocation of integer and real workspaces
      this->nprof = std::max(number_jacobian_nonzeros + this->number_variables, 5 /* heuristic */ * number_jacobian_nonzeros);
      this->mxws = this->compute_mxws();
      // 1 pointer hidden in lws (the quadratic program)
      constexpr size_t hidden_pointers_size = sizeof(intptr_t);
      this->mxlws = hidden_pointers_size + static_cast<size_t>(this->kmax) /* (required by bqpd.f) */ +
         9 * number_variables + number_constraints /* (required by sparseL.f) */;
      this->ws.resize(this->mxws);
      this->lws.resize(this->mxlws);
   }

   size_t BQPDSolver::compute_mxws() const {
      return static_cast<size_t>(this->kmax * (this->kmax + 9) / 2) + 2 * this->number_variables +
         this->number_constraints /* (required by bqpd.f) */ +
         5 * this->number_variables + this->nprof; /* (required by sparseL.f) */
   }

   QuadraticProgram& BQPDSolver::get_quadratic_program() {
      return *this->quadratic_program;
   }

   void BQPDSolver::solve(Statistics& /*statistics*/, const Vector<double>& initial_point, Direction& direction,
         const WarmstartInformation& warmstart_information) {
      // lazily allocate the BQPD scratch if the quadratic program was built directly from data (i.e.
      // initialize_memory(subproblem) was not called). The full solver path allocates it up front, so this
      // check is a no-op there.
      const size_t number_variables = this->quadratic_program->number_variables;
      const size_t number_constraints = this->quadratic_program->number_constraints;
      if (this->active_set.size() != number_variables + number_constraints) {
         this->allocate_workspace(number_variables, number_constraints, this->quadratic_program->number_jacobian_nonzeros,
            this->quadratic_program->has_curvature());
      }

      // initialize wsc_ common block (Hessian & workspace for BQPD)
      // setting the common block here ensures that several instances of BQPD can run simultaneously
      WSC.mxws = static_cast<int>(this->mxws);
      WSC.mxlws = static_cast<int>(this->mxlws);
      WSC.kk = 0; // length of ws that is used by gdotx
      WSC.ll = 0; // length of lws that is used by gdotx

      // hide a pointer to the quadratic program so that the gdotx callback can compute Hessian-vector products
      hide_pointer(0, this->lws.data(), *this->quadratic_program);
      WSC.ll += sizeof(intptr_t);

      if (this->print_subproblem) {
         this->display_subproblem(initial_point);
      }
      this->solve_subproblem(initial_point, direction, warmstart_information);
   }

   SolverWorkspace& BQPDSolver::get_workspace() {
      return this->quadratic_program->get_workspace();
   }

   // private member functions

   void BQPDSolver::display_subproblem(const Vector<double>& initial_point) const {
      const BQPDQuadraticProgram& quadratic_program = *this->quadratic_program;
      const size_t number_jacobian_nonzeros = quadratic_program.workspace.jacobian_values.size();
      DEBUG << "Subproblem:\n";
      DEBUG << "Linear objective part: " << view(quadratic_program.workspace.gradients, 0, quadratic_program.number_variables) << '\n';
      // note: Hessian values may not be available yet
      DEBUG << "Jacobian: " << view(quadratic_program.workspace.gradients, quadratic_program.number_variables,
         quadratic_program.number_variables + number_jacobian_nonzeros) << '\n';
      for (size_t variable_index: Range(quadratic_program.number_variables)) {
         DEBUG << "d" << variable_index << " in [" << quadratic_program.lower_bounds[variable_index] << ", " <<
            quadratic_program.upper_bounds[variable_index] << "]\n";
      }
      for (size_t constraint_index: Range(quadratic_program.number_constraints)) {
         DEBUG << "linearized c" << constraint_index << " in [" <<
            quadratic_program.lower_bounds[quadratic_program.number_variables + constraint_index] << ", " <<
            quadratic_program.upper_bounds[quadratic_program.number_variables + constraint_index] << "]\n";
      }
      DEBUG << "Initial point: " << initial_point << '\n';
   }

   // we use (x*4)/3 to increase by 33% TODO add an option
   bool BQPDSolver::check_sufficient_workspace_size(BQPDStatus bqpd_status) {
      switch (bqpd_status) {
         case BQPDStatus::REDUCED_HESSIAN_INSUFFICIENT_SPACE:
            // increase kmax
            this->kmax = (this->kmax*4)/3;
            return false;

         case BQPDStatus::SPARSE_INSUFFICIENT_SPACE:
            // allocate more size for (sparse) factors
            this->nprof *= 2;
            this->mxws = this->compute_mxws();
            // this->mxws = (this->mxws*4)/3;
            this->mxlws = (this->mxlws*4)/3;
            this->ws.resize(this->mxws);
            this->lws.resize(this->mxlws);
            WSC.mxws = static_cast<int>(this->mxws);
            WSC.mxlws = static_cast<int>(this->mxlws);
            return false;

         default:
            break;
      }
      return true;
   }

   void BQPDSolver::solve_subproblem(const Vector<double>& initial_point, Direction& direction,
         const WarmstartInformation& warmstart_information) {
      BQPDQuadraticProgram& quadratic_program = *this->quadratic_program;
      view(direction.primals, 0, quadratic_program.number_variables) = view(initial_point, 0, quadratic_program.number_variables);
      const int n = static_cast<int>(quadratic_program.number_variables);
      const int m = static_cast<int>(quadratic_program.number_constraints);

      const BQPDMode mode = BQPDSolver::determine_mode(warmstart_information);
      const int mode_integer = static_cast<int>(mode);

      // solve the LP/QP
      bool termination = false;
      while (!termination) {
         DEBUG2 << "Running BQPD\n";
         BQPD(&n, &m, &this->k, &this->kmax, quadratic_program.workspace.gradients.data(),
            quadratic_program.workspace.gradients_sparsity.data(), direction.primals.data(),
            quadratic_program.lower_bounds.data(), quadratic_program.upper_bounds.data(), &direction.subproblem_objective,
            &this->fmin, this->gradient_solution.data(), this->residuals.data(), this->w.data(), this->e.data(),
            this->active_set.data(), this->alp.data(), this->lp.data(), &this->mlp, &this->peq_solution,
            this->ws.data(), this->lws.data(), &mode_integer, &this->ifail, this->info.data(), &this->iprint, &this->nout);
         DEBUG2 << "Ran BQPD\n";
         const BQPDStatus bqpd_status = BQPDSolver::bqpd_status_from_int(this->ifail);
         termination = this->check_sufficient_workspace_size(bqpd_status);
         if (termination) {
            direction.status = BQPDSolver::status_from_bqpd_status(bqpd_status);
         }
      }

      // project primal solution into bounds
      for (size_t variable_index: Range(quadratic_program.number_variables)) {
         direction.primals[variable_index] = std::min(std::max(direction.primals[variable_index],
            quadratic_program.lower_bounds[variable_index]), quadratic_program.upper_bounds[variable_index]);
      }
      // gather the multipliers (the dual-displacement mapping is performed by IQPSolver)
      this->set_multipliers(quadratic_program.number_variables, direction.multipliers);
   }

   BQPDMode BQPDSolver::determine_mode(const WarmstartInformation& warmstart_information) {
      BQPDMode mode = BQPDMode::USER_DEFINED_ACTIVE_SET;
      // if problem structure changed, use cold start
      if (warmstart_information.hessian_sparsity_changed || warmstart_information.jacobian_sparsity_changed) {
         mode = BQPDMode::ACTIVE_SET_EQUALITIES;
      }
      // if only the variable bounds changed, reuse the active set estimate and the Jacobian information
      else if (warmstart_information.trust_region_changed && !warmstart_information.new_iterate &&
            !warmstart_information.constraint_bounds_changed) {
         mode = BQPDMode::UNCHANGED_ACTIVE_SET_AND_JACOBIAN;
      }
      return mode;
   }

   void BQPDSolver::set_multipliers(size_t number_variables, Multipliers& direction_multipliers) const {
      const BQPDQuadraticProgram& quadratic_program = *this->quadratic_program;
      direction_multipliers.reset();
      // active constraints
      for (size_t active_constraint_index: Range(number_variables - static_cast<size_t>(this->k))) {
         const size_t index = static_cast<size_t>(std::abs(this->active_set[active_constraint_index])) - Indexing::Fortran_indexing;

         // bound constraint
         if (index < number_variables) {
            // variable is not fixed
            if (quadratic_program.lower_bounds[index] < quadratic_program.upper_bounds[index]) {
               if (0 <= this->active_set[active_constraint_index]) { // lower bound active
                  direction_multipliers.lower_bounds[index] = this->residuals[index];
               }
               else { // upper bound active
                  direction_multipliers.upper_bounds[index] = -this->residuals[index];
               }
            }
            else { // variable fixed
               if (0. <= this->residuals[index]) {
                  direction_multipliers.lower_bounds[index] = this->residuals[index];
               }
               else {
                  direction_multipliers.upper_bounds[index] = this->residuals[index];
               }
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
            DEBUG << "BQPD error: insufficient LP space\n";
            return SubproblemStatus::ERROR;
         case BQPDStatus::MAX_RESTARTS_REACHED:
            DEBUG << "BQPD max restarts reached\n";
            return SubproblemStatus::ERROR;
         default:
            throw std::invalid_argument("The BQPD ifail is not consistent with the Uno status values");
      }
   }
} // namespace

void hessian_vector_product(int* dimension, const double vector[], const double /*ws*/[], const int lws[], double result[]) {
   assert(dimension != nullptr && "BQPDSolver::hessian_vector_product: the dimension n passed by pointer is NULL");

   // retrieve the quadratic program and delegate the (Subproblem-free) Hessian-vector product
   uno::BQPDQuadraticProgram* quadratic_program = uno::retrieve_pointer<uno::BQPDQuadraticProgram>(0, lws);
   assert(quadratic_program != nullptr);
   quadratic_program->compute_hessian_vector_product(*dimension, vector, result);
}
