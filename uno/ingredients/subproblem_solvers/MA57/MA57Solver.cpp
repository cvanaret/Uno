// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cassert>
#include <optional>
#include <utility>
#include <vector>
#include "MA57Solver.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Direction.hpp"
#include "tools/Logger.hpp"
#include "fortran_interface.h"

#define MA57_set_default_parameters FC_GLOBAL(ma57id, MA57ID)
#define MA57_symbolic_analysis FC_GLOBAL(ma57ad, MA57AD)
#define MA57_numerical_factorization FC_GLOBAL(ma57bd, MA57BD)
#define MA57_linear_solve FC_GLOBAL(ma57cd, MA57CD)
#define MA57_linear_solve_with_iterative_refinement FC_GLOBAL(ma57dd, MA57DD)
#define MA57_enlarge_workspace FC_GLOBAL(ma57ed, MA57ED)

namespace uno {
   extern "C" {
      void MA57_set_default_parameters(double cntl[], int icntl[]);

      void MA57_symbolic_analysis(const int* n, const int* ne, const int irn[], const int jcn[], const int* lkeep,
         /*const*/ int keep[], int iwork[], /* const */ int icntl[], int info[], double rinfo[]);

      void MA57_numerical_factorization(const int* n, int* ne, const double a[], double fact[], const int* lfact,
         int ifact[], const int* lifact, const int* lkeep, const int keep[], int iwork[], int icntl[], double cntl[],
         int info[], double rinfo[]);

      void MA57_linear_solve(const int* job, const int* n, double fact[], int* lfact, int ifact[], int* lifact, const int* nrhs,
         double rhs[], const int* lrhs, double work[], int* lwork, int iwork[], int icntl[], int info[]);

      void MA57_linear_solve_with_iterative_refinement(const int* job, const int* n, int* ne, const double a[], const int irn[],
         const int jcn[], double fact[], int* lfact, int ifact[], int* lifact, const double rhs[], double x[], double resid[],
         double work[], int iwork[], int icntl[], double cntl[], int info[], double rinfo[]);

      // enlarging of workspaces when numerical factorization runs out of memory
      void MA57_enlarge_workspace(const int* n, const int* ic, int keep[], const double fact[], const int* lfact,
         double newfac[], const int* lnew, const int ifact[], const int* lifact, int newifc[], const int* linew,
         int info[]);
   }

   namespace {
      bool is_error_code_insufficient_real_workspace(int error_code) {
         return error_code == 10 || error_code == -3;
      }

      bool is_error_code_insufficient_integer_workspace(int error_code) {
         return error_code == 11 || error_code == -4;
      }

      int get_larger_workspace_size(int size_current, std::optional<int> size_estimate) {
         const int size_new = std::max(size_current + 1, size_estimate.value_or(size_current));
         const int oversize_denominator_with_estimate = 5;  // add 20% on top of the estimate
         const int oversize_denominator_without_estimate = 2;  // grow by 50% if we don't have an estimate
         const int oversize_denominator = size_estimate.has_value() ? oversize_denominator_with_estimate : oversize_denominator_without_estimate;
         return size_new + size_new / oversize_denominator;
      }

      int get_larger_real_workspace_size(const MA57Workspace& workspace) {
         const bool have_estimate = (workspace.info[0] == -3);
         const auto lfact_estimate = have_estimate ? std::optional<int>{workspace.info[16]} : std::nullopt;
         return get_larger_workspace_size(workspace.lfact, lfact_estimate);
      }

      int get_larger_integer_workspace_size(const MA57Workspace& workspace) {
         const bool have_estimate = (workspace.info[0] == -4);
         const auto lifact_estimate = have_estimate ? std::optional<int>{workspace.info[17]} : std::nullopt;
         return get_larger_workspace_size(workspace.lifact, lifact_estimate);
      }
   }  // anonymous namespace

   MA57Solver::MA57Solver(): DirectSymmetricIndefiniteLinearSolver() {
      // set the default values of the controlling parameters
      MA57_set_default_parameters(this->workspace.cntl.data(), this->workspace.icntl.data());
      // suppress warning messages
      this->workspace.icntl[4] = 0;
      // iterative refinement enabled
      this->workspace.icntl[8] = 1;
   }

   void MA57Solver::initialize_hessian(const Subproblem& subproblem) {
      this->evaluation_space.initialize_hessian(subproblem);

      // workspace
      const size_t dimension = subproblem.number_variables;
      this->workspace.n = static_cast<int>(dimension);
      this->workspace.nnz = static_cast<int>(this->evaluation_space.number_matrix_nonzeros);
      this->workspace.lkeep = static_cast<int>(5 * dimension + this->evaluation_space.number_matrix_nonzeros +
         std::max(dimension, this->evaluation_space.number_matrix_nonzeros) + 42);
      this->workspace.keep.resize(static_cast<size_t>(this->workspace.lkeep));
      this->workspace.iwork.resize(5 * dimension);
      this->workspace.lwork = static_cast<int>(1.2 * static_cast<double>(dimension));
      this->workspace.work.resize(static_cast<size_t>(this->workspace.lwork));
      this->workspace.residuals.resize(dimension);
   }

   void MA57Solver::initialize_augmented_system(const Subproblem& subproblem) {
      this->evaluation_space.initialize_augmented_system(subproblem);

      // workspace
      const size_t dimension = subproblem.number_variables + subproblem.number_constraints;
      this->workspace.n = static_cast<int>(dimension);
      this->workspace.nnz = static_cast<int>(this->evaluation_space.number_matrix_nonzeros);
      this->workspace.lkeep = static_cast<int>(5 * dimension + this->evaluation_space.number_matrix_nonzeros +
         std::max(dimension, this->evaluation_space.number_matrix_nonzeros) + 42);
      this->workspace.keep.resize(static_cast<size_t>(this->workspace.lkeep));
      this->workspace.iwork.resize(5 * dimension);
      this->workspace.lwork = static_cast<int>(1.2 * static_cast<double>(dimension));
      this->workspace.work.resize(static_cast<size_t>(this->workspace.lwork));
      this->workspace.residuals.resize(dimension);
   }

   void MA57Solver::do_symbolic_analysis() {
      assert(!this->analysis_performed);

      // symbolic analysis
      MA57_symbolic_analysis(&this->workspace.n, &this->workspace.nnz, this->evaluation_space.matrix_row_indices.data(),
         this->evaluation_space.matrix_column_indices.data(), &this->workspace.lkeep, this->workspace.keep.data(),
         this->workspace.iwork.data(), this->workspace.icntl.data(), this->workspace.info.data(), this->workspace.rinfo.data());

      assert(0 <= this->workspace.info[0] && "MA57: the symbolic analysis failed");
      if (0 < this->workspace.info[0]) {
         WARNING << "MA57 has issued a warning: info(1) = " << workspace.info[0] << '\n';
      }

      // get LFACT and LIFACT and resize FACT and IFACT (no effect if resized to <= size)
      int lfact = 2 * this->workspace.info[8];
      int lifact = 2 * this->workspace.info[9];
      this->workspace.fact.resize(static_cast<size_t>(lfact));
      this->workspace.ifact.resize(static_cast<size_t>(lifact));
      // store the sizes of the symbolic analysis
      this->workspace.lfact = lfact;
      this->workspace.lifact = lifact;
      this->analysis_performed = true;
   }

   void MA57Solver::do_numerical_factorization(const double* matrix_values) {
      assert(this->analysis_performed);

      bool factorization_done = false;
      while (!factorization_done) {
         // numerical factorization
         MA57_numerical_factorization(&this->workspace.n, &this->workspace.nnz, matrix_values, this->workspace.fact.data(),
            &this->workspace.lfact, this->workspace.ifact.data(), &this->workspace.lifact, &this->workspace.lkeep,
            this->workspace.keep.data(), this->workspace.iwork.data(), this->workspace.icntl.data(), this->workspace.cntl.data(),
            this->workspace.info.data(), this->workspace.rinfo.data());

         if (is_error_code_insufficient_real_workspace(this->workspace.info[0]) ||
             is_error_code_insufficient_integer_workspace(this->workspace.info[0])) {
            const bool is_real_workspace = is_error_code_insufficient_real_workspace(this->workspace.info[0]);

            const int lnewfact = !is_real_workspace ? 0 : get_larger_real_workspace_size(this->workspace);
            const int lnewifact = is_real_workspace ? 0 : get_larger_integer_workspace_size(this->workspace);
            std::vector<double> newfact(static_cast<size_t>(lnewfact));
            std::vector<int> newifact(static_cast<size_t>(lnewifact));
            const int enlarge_target = is_real_workspace ? 0 : 1;

            MA57_enlarge_workspace(&this->workspace.n, &enlarge_target, this->workspace.keep.data(), this->workspace.fact.data(),
               &this->workspace.lfact, newfact.data(), &lnewfact, this->workspace.ifact.data(), &this->workspace.lifact,
               newifact.data(), &lnewifact, this->workspace.info.data());

            if (is_real_workspace) {
               this->workspace.fact = std::move(newfact);
               this->workspace.lfact = lnewfact;
            } else {
               this->workspace.ifact = std::move(newifact);
               this->workspace.lifact = lnewifact;
            }
         } else {
            factorization_done = true;
         }
      }
      this->factorization_performed = true;
   }

   void MA57Solver::solve_indefinite_system(const Vector<double>& matrix_values, const Vector<double>& rhs, Vector<double>& result) {
      assert(this->factorization_performed);

      // solve
      const int lrhs = this->workspace.n; // integer, length of rhs

      // solve the linear system
      if (this->use_iterative_refinement) {
         MA57_linear_solve_with_iterative_refinement(&this->workspace.job, &this->workspace.n, &this->workspace.nnz,
            matrix_values.data(), this->evaluation_space.matrix_row_indices.data(), this->evaluation_space.matrix_column_indices.data(),
            this->workspace.fact.data(), &this->workspace.lfact, this->workspace.ifact.data(), &this->workspace.lifact,
            rhs.data(), result.data(), this->workspace.residuals.data(), this->workspace.work.data(), this->workspace.iwork.data(),
            this->workspace.icntl.data(), this->workspace.cntl.data(), this->workspace.info.data(), this->workspace.rinfo.data());
      }
      else {
         // copy rhs into result (overwritten by MA57)
         result = rhs;

         MA57_linear_solve(&this->workspace.job, &this->workspace.n, this->workspace.fact.data(), &this->workspace.lfact,
            this->workspace.ifact.data(), &this->workspace.lifact, &this->workspace.nrhs, result.data(), &lrhs,
            this->workspace.work.data(), &this->workspace.lwork, this->workspace.iwork.data(), this->workspace.icntl.data(),
            this->workspace.info.data());
      }
   }

   void MA57Solver::solve_indefinite_system(Statistics& statistics, const Subproblem& subproblem, Direction& direction,
         const WarmstartInformation& warmstart_information) {
      // set up the linear system by evaluating the functions at the current iterate
      this->evaluation_space.set_up_linear_system(statistics, subproblem, *this, warmstart_information);
      // solve the linear system
      this->solve_indefinite_system(this->evaluation_space.matrix_values, this->evaluation_space.rhs, this->evaluation_space.solution);
      // assemble the full primal-dual direction
      subproblem.assemble_primal_dual_direction(this->evaluation_space.solution, direction);
      if (this->matrix_is_singular()) {
         direction.status = SubproblemStatus::INFEASIBLE;
      }
   }

   Inertia MA57Solver::get_inertia() const {
      // rank = number_positive_eigenvalues + number_negative_eigenvalues
      // n = rank + number_zero_eigenvalues
      const size_t rank = this->rank();
      const size_t number_negative_eigenvalues = this->number_negative_eigenvalues();
      const size_t number_positive_eigenvalues = rank - number_negative_eigenvalues;
      const size_t number_zero_eigenvalues = static_cast<size_t>(this->workspace.n) - rank;
      return {number_positive_eigenvalues, number_negative_eigenvalues, number_zero_eigenvalues};
   }

   size_t MA57Solver::number_negative_eigenvalues() const {
      return static_cast<size_t>(this->workspace.info[23]);
   }

   /*
   bool MA57Solver::matrix_is_positive_definite() const {
      // positive definite = non-singular and no negative eigenvalues
      return not this->matrix_is_singular() && this->number_negative_eigenvalues() == 0;
   }
   */

   bool MA57Solver::matrix_is_singular() const {
      return (this->workspace.info[0] == 4);
   }

   size_t MA57Solver::rank() const {
      return static_cast<size_t>(this->workspace.info[24]);
   }

   EvaluationSpace& MA57Solver::get_evaluation_space() {
      return this->evaluation_space;
   }
} // namespace