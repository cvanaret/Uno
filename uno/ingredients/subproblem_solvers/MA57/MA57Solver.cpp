// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cassert>
#include <optional>
#include <utility>
#include <vector>
#include "MA57Solver.hpp"
#include "tools/Logger.hpp"

#ifdef HSL_RUNTIME_LOADING
#include <stdexcept>
#include "ingredients/subproblem_solvers/HSL/HSLLoader.hpp"
// route the calls through the runtime-resolved function pointers
#define MA57_set_default_parameters uno::hsl_ma57id
#define MA57_symbolic_analysis uno::hsl_ma57ad
#define MA57_numerical_factorization uno::hsl_ma57bd
#define MA57_linear_solve uno::hsl_ma57cd
#define MA57_enlarge_workspace uno::hsl_ma57ed
#else
#include "fortran_interface.h"
#define MA57_set_default_parameters FC_GLOBAL(ma57id, MA57ID)
#define MA57_symbolic_analysis FC_GLOBAL(ma57ad, MA57AD)
#define MA57_numerical_factorization FC_GLOBAL(ma57bd, MA57BD)
#define MA57_linear_solve FC_GLOBAL(ma57cd, MA57CD)
#define MA57_enlarge_workspace FC_GLOBAL(ma57ed, MA57ED)
#endif
#if defined(HAS_HSL) && !defined(HSL_RUNTIME_LOADING)
// libhsl.h declares some functions with a parameter called "new", which is a reserved C++ keyword.
// Temporarily rename it with #define
#define new hsl_new
#include <libhsl.h>
#undef new
#endif

namespace uno {
#ifndef HSL_RUNTIME_LOADING
   extern "C" {
      void MA57_set_default_parameters(double cntl[], int icntl[]);

      void MA57_symbolic_analysis(const int* n, const int* ne, const int irn[], const int jcn[], const int* lkeep,
         /*const*/ int keep[], int iwork[], /* const */ int icntl[], int info[], double rinfo[]);

      void MA57_numerical_factorization(const int* n, int* ne, const double a[], double fact[], const int* lfact,
         int ifact[], const int* lifact, const int* lkeep, const int keep[], int iwork[], int icntl[], double cntl[],
         int info[], double rinfo[]);

      void MA57_linear_solve(const int* job, const int* n, double fact[], int* lfact, int ifact[], int* lifact, const int* nrhs,
         double rhs[], const int* lrhs, double work[], int* lwork, int iwork[], int icntl[], int info[]);

      // enlarging of workspaces when numerical factorization runs out of memory
      void MA57_enlarge_workspace(const int* n, const int* ic, int keep[], const double fact[], const int* lfact,
         double newfac[], const int* lnew, const int ifact[], const int* lifact, int newifc[], const int* linew,
         int info[]);
   }
#endif

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
#ifdef HSL_RUNTIME_LOADING
      if (!ma57_symbols_available()) {
         throw std::runtime_error("Uno: the MA57 solver was requested but the HSL library could not be loaded at "
            "runtime (set the UNO_HSL_LIBRARY environment variable to point at libhsl)");
      }
#endif
#if defined(HAS_HSL) && !defined(HSL_RUNTIME_LOADING)
      INFO << "Running MA57 v" << LIBHSL_VER_MAJOR << "." << LIBHSL_VER_MINOR << "." << LIBHSL_VER_PATCH << '\n';
#else
      INFO << "Running MA57 v1.0.0\n";
#endif
      // set the default values of the controlling parameters
      MA57_set_default_parameters(this->workspace.cntl.data(), this->workspace.icntl.data());
      MA57_ICNTL(5) = 0; // suppress warning messages
      MA57_ICNTL(6) = 5; // pivot order (auto between METIS and MC47)
      MA57_ICNTL(7) = 1; // numerical pivoting (uses threshold in CNTL(1))
      MA57_ICNTL(11) = 16; // block size used by the Level 3 BLAS
      MA57_ICNTL(12) = 16; // assembly tree
      MA57_ICNTL(15) = MA57Settings::mc64_scaling; // MC64 scaling (disabled)
      MA57_ICNTL(16) = 0; // small entries removed (disabled)
      MA57_CNTL(1) = MA57Settings::pivoting_threshold; // pivoting threshold
   }

   void MA57Solver::initialize_memory() {
      this->workspace.n = static_cast<int>(this->linear_system.dimension);
      this->workspace.nnz = static_cast<int>(this->linear_system.number_nonzeros);
      this->workspace.lkeep = static_cast<int>(5 * this->linear_system.dimension + this->linear_system.number_nonzeros +
         std::max(this->linear_system.dimension, this->linear_system.number_nonzeros) + 42);
      this->workspace.keep.resize(static_cast<size_t>(this->workspace.lkeep));
      this->workspace.iwork.resize(5 * this->linear_system.dimension);
      this->workspace.lwork = static_cast<int>(static_cast<double>(this->linear_system.dimension));
      this->workspace.work.resize(static_cast<size_t>(this->workspace.lwork));
   }

   void MA57Solver::do_symbolic_analysis() {
      assert(!this->analysis_performed);

      MA57_symbolic_analysis(&this->workspace.n, &this->workspace.nnz, this->linear_system.matrix_row_indices.data(),
         this->linear_system.matrix_column_indices.data(), &this->workspace.lkeep, this->workspace.keep.data(),
         this->workspace.iwork.data(), this->workspace.icntl.data(), this->workspace.info.data(), this->workspace.rinfo.data());

      if (MA57_INFO(1) < 0) {
         throw std::runtime_error("MA57: the symbolic analysis failed");
      }
      if (0 < MA57_INFO(1)) {
         WARNING << "MA57 has issued a warning: info(1) = " << MA57_INFO(1) << '\n';
      }

      // get LFACT and LIFACT and resize FACT and IFACT (no effect if resized to <= size)
      this->workspace.lfact = static_cast<int>(MA57Settings::allocation_safety_factor * MA57_INFO(11));
      this->workspace.lifact = static_cast<int>(MA57Settings::allocation_safety_factor * MA57_INFO(12));
      this->workspace.fact.resize(static_cast<size_t>(this->workspace.lfact));
      this->workspace.ifact.resize(static_cast<size_t>(this->workspace.lifact));
      this->analysis_performed = true;
   }

   void MA57Solver::do_numerical_factorization(bool /*is_matrix_positive_definite*/) {
      assert(this->analysis_performed);

      bool factorization_done = false;
      while (!factorization_done) {
         // numerical factorization
         MA57_numerical_factorization(&this->workspace.n, &this->workspace.nnz, this->linear_system.matrix_values.data(),
            this->workspace.fact.data(), &this->workspace.lfact, this->workspace.ifact.data(), &this->workspace.lifact,
            &this->workspace.lkeep, this->workspace.keep.data(), this->workspace.iwork.data(), this->workspace.icntl.data(),
            this->workspace.cntl.data(), this->workspace.info.data(), this->workspace.rinfo.data());

         if (is_error_code_insufficient_real_workspace(MA57_INFO(1)) || is_error_code_insufficient_integer_workspace(MA57_INFO(1))) {
            const bool is_real_workspace = is_error_code_insufficient_real_workspace(this->workspace.info[0]);

            const int lnewfact = !is_real_workspace ? 0 : get_larger_real_workspace_size(this->workspace);
            const int lnewifact = is_real_workspace ? 0 : get_larger_integer_workspace_size(this->workspace);
            std::vector<double> newfact(static_cast<size_t>(lnewfact));
            std::vector<int> newifact(static_cast<size_t>(lnewifact));
            const int enlarge_target = is_real_workspace ? 0 : 1;

            DEBUG << "Enlarging the MA57 workspace\n";
            MA57_enlarge_workspace(&this->workspace.n, &enlarge_target, this->workspace.keep.data(), this->workspace.fact.data(),
               &this->workspace.lfact, newfact.data(), &lnewfact, this->workspace.ifact.data(), &this->workspace.lifact,
               newifact.data(), &lnewifact, this->workspace.info.data());

            if (is_real_workspace) {
               this->workspace.fact = std::move(newfact);
               this->workspace.lfact = lnewfact;
            }
            else {
               this->workspace.ifact = std::move(newifact);
               this->workspace.lifact = lnewifact;
            }
         }
         else if (MA57_INFO(1) < 0) {
            throw std::runtime_error("MA57 fatal error");
         }
         else {
            factorization_done = true;
         }
      }
      this->factorization_performed = true;
   }

   void MA57Solver::solve_indefinite_system(double* solution) {
      return this->solve_indefinite_system(this->linear_system.rhs.data(), solution, 1);
   }

   void MA57Solver::solve_indefinite_system(const double* rhs, double* solution, size_t number_of_rhs) {
      assert(this->factorization_performed);

      const int nrhs = static_cast<int>(number_of_rhs);
      const int lrhs = this->workspace.n; // leading dimension of the right-hand side block

      // copy the rhs block into the solution block (overwritten by MA57)
      const size_t dimension = static_cast<size_t>(this->workspace.n);
      view(solution, number_of_rhs * dimension) = view(rhs, number_of_rhs * dimension);

      // ma57cd requires a work array of size n*nrhs; enlarge it if needed
      const int required_lwork = this->workspace.n * nrhs;
      if (this->workspace.lwork < required_lwork) {
         this->workspace.lwork = required_lwork;
         this->workspace.work.resize(static_cast<size_t>(this->workspace.lwork));
      }

      MA57_linear_solve(&this->workspace.job, &this->workspace.n, this->workspace.fact.data(), &this->workspace.lfact,
         this->workspace.ifact.data(), &this->workspace.lifact, &nrhs, solution, &lrhs,
         this->workspace.work.data(), &this->workspace.lwork, this->workspace.iwork.data(), this->workspace.icntl.data(),
         this->workspace.info.data());
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
      return static_cast<size_t>(MA57_INFO(24));
   }

   /*
   bool MA57Solver::matrix_is_positive_definite() const {
      // positive definite = non-singular and no negative eigenvalues
      return not this->matrix_is_singular() && this->number_negative_eigenvalues() == 0;
   }
   */

   bool MA57Solver::matrix_is_singular() const {
      return (MA57_INFO(1) == 4);
   }

   size_t MA57Solver::rank() const {
      return static_cast<size_t>(this->workspace.info[24]);
   }

   LinearSystem& MA57Solver::get_linear_system() {
      return this->linear_system;
   }

   COOLinearSystem& MA57Solver::get_coo_linear_system() {
      return this->linear_system;
   }

   // protected member functions

   int& MA57Solver::MA57_ICNTL(size_t index) {
      // handle the Fortran indexing (starting at 1)
      return this->workspace.icntl[index-1];
   }

   double& MA57Solver::MA57_CNTL(size_t index) {
      // handle the Fortran indexing (starting at 1)
      return this->workspace.cntl[index-1];
   }

   int MA57Solver::MA57_INFO(size_t index) const {
      // handle the Fortran indexing (starting at 1)
      return this->workspace.info[index-1];
   }
} // namespace