// Copyright (c) 2026 Uno contributors
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include "BQPDCPPSolver.hpp"
#include "ingredients/subproblem_solvers/BQPD/BQPDQuadraticProgram.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/Indexing.hpp"
#include "linear_algebra/Vector.hpp"
#include "linear_algebra/View.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Multipliers.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "symbolic/Range.hpp"
#include "tools/Logger.hpp"
// BQPD-CPP (the C++ port) public API
#include "bqpd/BqpdSolver.hpp"
#include "bqpd/Common.hpp"
#include "bqpd/HessianProduct.hpp"
#include "bqpd/SparseA.hpp"
#include "bqpd/SparseFactor.hpp"

namespace uno {
   namespace {
      // Adapter: BQPD-CPP asks for v := G*x through HessianProduct::hprod; delegate to the
      // reused BQPDQuadraticProgram (same symmetric COO matvec the Fortran gdotx used).
      struct HessianProductAdapter : public bqpd::HessianProduct {
         const BQPDQuadraticProgram* quadratic_program;
         bool zero;
         HessianProductAdapter(const BQPDQuadraticProgram* qp, bool is_zero): quadratic_program(qp), zero(is_zero) {}
         void hprod(int dimension, const bqpd::real* vector, bqpd::real* result) const override {
            this->quadratic_program->compute_hessian_vector_product(dimension, vector, result);
         }
         bool is_zero() const override { return this->zero; }
      };

      // Same logic as BQPDSolver::determine_mode (bqpd modes 0..6).
      int determine_mode(const WarmstartInformation& warmstart_information) {
         int mode = 2; // USER_DEFINED_ACTIVE_SET (hot start)
         if (warmstart_information.hessian_sparsity_changed || warmstart_information.jacobian_sparsity_changed) {
            mode = 1; // ACTIVE_SET_EQUALITIES (cold start)
         }
         else if (warmstart_information.trust_region_changed && !warmstart_information.new_iterate &&
               !warmstart_information.constraint_bounds_changed) {
            mode = 4; // UNCHANGED_ACTIVE_SET_AND_JACOBIAN (BQPD-CPP clamps this to 2)
         }
         return mode;
      }

      // Map BQPD-CPP's Ifail (integer values identical to Fortran ifail 0..8) to a SubproblemStatus.
      SubproblemStatus status_from_ifail(int ifail) {
         switch (ifail) {
            case 0: return SubproblemStatus::OPTIMAL;             // KT point
            case 1: return SubproblemStatus::UNBOUNDED_PROBLEM;
            case 3: return SubproblemStatus::INFEASIBLE;
            case 2:                                               // bound inconsistency
            case 4:                                               // incorrect parameter
            case 5:                                               // insufficient LP space
            case 8:                                               // max restarts reached
               DEBUG << "BQPD-CPP error/limit: ifail = " << ifail << '\n';
               return SubproblemStatus::ERROR;
            default:
               throw std::invalid_argument("The BQPD-CPP ifail is not consistent with the Uno status values");
         }
      }
   } // namespace

   BQPDCPPSolver::BQPDCPPSolver(const Options& options):
         QPSolver(),
         alp(static_cast<std::size_t>(this->mlp)),
         lp(static_cast<std::size_t>(this->mlp)),
         print_subproblem(options.get_bool("print_subproblem")) {
      this->quadratic_program = std::make_unique<BQPDQuadraticProgram>();
      this->ctx.prt.iprint = 0; // silence the port's trace; UNO drives its own logging
      this->ctx.eps.emin = 0.0;
   }

   BQPDCPPSolver::~BQPDCPPSolver() = default;

   void BQPDCPPSolver::initialize_memory(const Subproblem& subproblem) {
      this->quadratic_program = std::make_unique<BQPDQuadraticProgram>();
      this->quadratic_program->initialize_memory(subproblem);
      this->allocate_workspace(subproblem.number_variables, subproblem.number_constraints,
         subproblem.number_jacobian_nonzeros(), subproblem.has_curvature());
   }

   void BQPDCPPSolver::allocate_workspace(std::size_t number_variables, std::size_t number_constraints,
         std::size_t number_jacobian_nonzeros, bool is_qp) {
      const std::size_t nm = number_variables + number_constraints;
      this->w.resize(nm);
      this->gradient_solution.resize(number_variables);
      this->residuals.resize(nm);
      this->e.resize(nm);

      // default active set: every variable/constraint index (Fortran 1-based)
      this->active_set.resize(nm);
      for (std::size_t index: Range(nm)) {
         this->active_set[index] = static_cast<int>(index + Indexing::Fortran_indexing);
      }

      // kmax (max nullspace size) is 0 for an LP
      this->kmax = is_qp ? std::min(static_cast<int>(number_variables), 500) : 0;
      // nprof: the sparse row-spike profile; unknown a priori, grows on ifail = 7
      this->nprof = std::max(number_jacobian_nonzeros + number_variables, 5 * number_jacobian_nonzeros);
   }

   QuadraticProgram& BQPDCPPSolver::get_quadratic_program() {
      return *this->quadratic_program;
   }

   SolverWorkspace& BQPDCPPSolver::get_workspace() {
      return *this->quadratic_program;
   }

   void BQPDCPPSolver::solve(Statistics& /*statistics*/, const Vector<double>& initial_point, Direction& direction,
         const WarmstartInformation& warmstart_information) {
      BQPDQuadraticProgram& qp = *this->quadratic_program;

      // lazily allocate scratch if the QP was built directly from data (no initialize_memory)
      if (this->active_set.size() != qp.number_variables + qp.number_constraints) {
         this->allocate_workspace(qp.number_variables, qp.number_constraints,
            qp.number_jacobian_nonzeros, qp.has_curvature());
      }

      const int n = static_cast<int>(qp.number_variables);
      const int m = static_cast<int>(qp.number_constraints);
      const int nm = n + m;

      // BQPD-CPP consumes exactly the Fortran la/a layout BQPDQuadraticProgram already builds.
      std::vector<bqpd::real> a(qp.gradients.data(), qp.gradients.data() + qp.gradients.size());
      std::vector<int> la(qp.gradients_sparsity.data(), qp.gradients_sparsity.data() + qp.gradients_sparsity.size());

      // mode from warmstart info, then clamp >2 -> 2 (BQPD-CPP does not implement modes 3-6)
      int mode = determine_mode(warmstart_information);
      if (mode > 2) {
         mode = 2;
      }

      const bool is_qp = qp.has_curvature();

      // warm-start entry state to restore on each ifail=6/7 retry
      const int k_entry = this->k;
      const std::vector<int> active_set_entry = this->active_set;

      bqpd::Result result;
      bool termination = false;
      while (!termination) {
         this->k = k_entry;
         this->active_set = active_set_entry;

         // workspace sizing (mirrors qp_driver.cpp sparse path): bqpd needs kkk reals + lll ints,
         // sparseL.f needs 5n + nprof reals and 9n + m ints; nprof grows on ifail = 7.
         const int kkk = this->kmax * (this->kmax + 9) / 2 + nm + n;
         const int lll = this->kmax;
         const int mxws = kkk + 5 * n + static_cast<int>(this->nprof);
         const int mxlws = lll + 9 * n + m + nm;
         if (mxws != this->current_mxws || mxlws != this->current_mxlws) {
            this->workspace.resize(mxws, mxlws);
            this->current_mxws = mxws;
            this->current_mxlws = mxlws;
         }

         // (re)construct the seams over the shared workspace/context
         bqpd::SparseA A(n, m, a, la);
         bqpd::SparseFactor factor(A, this->workspace, this->ctx);
         HessianProductAdapter hessian(this->quadratic_program.get(), !is_qp);
         bqpd::BqpdSolver solver(A, factor, hessian, this->workspace, this->ctx);

         // x starts at the initial point
         view(direction.primals, 0, qp.number_variables) = view(initial_point, 0, qp.number_variables);

         result = solver.solve(n, m, this->k, this->kmax, this->mlp, static_cast<bqpd::Mode>(mode),
            direction.primals.data(), qp.lower_bounds.data(), qp.upper_bounds.data(), this->fmin,
            this->gradient_solution.data(), this->residuals.data(), this->w.data(), this->e.data(),
            this->active_set.data(), this->alp.data(), this->lp.data());

         const int ifail = static_cast<int>(result.ifail);
         if (ifail == 6) { // reduced Hessian too small: grow kmax and retry
            this->kmax = (this->kmax * 4) / 3;
         }
         else if (ifail == 7) { // sparse factors too small: grow nprof and retry
            this->nprof *= 2;
         }
         else {
            termination = true;
            direction.status = status_from_ifail(ifail);
         }
      }

      direction.subproblem_objective = result.f;
      this->peq = result.peq;

      // project primal solution into bounds
      for (std::size_t variable_index: Range(qp.number_variables)) {
         direction.primals[variable_index] = std::min(std::max(direction.primals[variable_index],
            qp.lower_bounds[variable_index]), qp.upper_bounds[variable_index]);
      }
      // gather the multipliers (the dual-displacement mapping is performed by IQPSolver)
      this->set_multipliers(qp.number_variables, direction.multipliers);
   }

   // identical to BQPDSolver::set_multipliers: multipliers live in residuals[abs(ls[i]) - 1]
   void BQPDCPPSolver::set_multipliers(std::size_t number_variables, Multipliers& direction_multipliers) const {
      const BQPDQuadraticProgram& qp = *this->quadratic_program;
      direction_multipliers.reset();
      for (std::size_t active_constraint_index: Range(number_variables - static_cast<std::size_t>(this->k))) {
         const std::size_t index = static_cast<std::size_t>(std::abs(this->active_set[active_constraint_index])) -
            Indexing::Fortran_indexing;

         // bound constraint
         if (index < number_variables) {
            if (qp.lower_bounds[index] < qp.upper_bounds[index]) { // variable not fixed
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
            const std::size_t constraint_index = index - number_variables;
            if (0 <= this->active_set[active_constraint_index]) { // lower bound active
               direction_multipliers.constraints[constraint_index] = this->residuals[index];
            }
            else { // upper bound active
               direction_multipliers.constraints[constraint_index] = -this->residuals[index];
            }
         }
      }
   }
} // namespace
