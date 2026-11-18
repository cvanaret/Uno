// Copyright (c) 2024-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cassert>
#include <limits>
#include <cmath>
#include <stdexcept>
#include "MUMPSSolver.hpp"
#include "linear_algebra/Vector.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#if defined(HAS_MPI) && defined(MUMPS_PARALLEL)
#include "mpi.h"
#endif

#define USE_COMM_WORLD (-987654)

namespace uno {
   MUMPSSolver::MUMPSSolver(bool use_scaling): DirectSymmetricIndefiniteLinearSolver(), use_scaling(use_scaling) {
      INFO << "Running MUMPS v" << MUMPS_VERSION << '\n';
      this->workspace.sym = MUMPSSolver::GENERAL_SYMMETRIC;
#if defined(HAS_MPI) && defined(MUMPS_PARALLEL)
      // TODO load number of processes from option file
      this->workspace.par = 1;
#else
      this->workspace.par = 1;
#endif
      this->workspace.job = MUMPSSolver::JOB_INIT;
      this->workspace.comm_fortran = USE_COMM_WORLD;
      dmumps_c(&this->workspace);
      // control parameters
      MUMPS_ICNTL(1) = -1; // output stream for error messages (off)
      MUMPS_ICNTL(2) = -1; // output stream for diagnostic printing (off)
      MUMPS_ICNTL(3) = -1; // output stream for global information (off)
      MUMPS_ICNTL(4) = 0; // level of printing for error, warning, and diagnostic message (off)

      MUMPS_ICNTL(6) = MUMPSSettings::permuting_scaling; // column permutation (none)
      MUMPS_ICNTL(7) = MUMPSSettings::pivot_order; // pivot order (AMD)
      MUMPS_ICNTL(8) = use_scaling ? MUMPSSettings::scaling : 0; // scaling strategy (diagonal scaling during factorization, or none)
      MUMPS_ICNTL(10) = 0; // iterative refinement (none)
      MUMPS_ICNTL(13) = 1; // parallelism of the root nod (ScaLAPACK not used, partly recover parallelism of the root node)
      MUMPS_ICNTL(14) = MUMPSSettings::mem_percent; // percentage increase in the estimated working space (35)
      MUMPS_ICNTL(24) = 1; // controls the detection of “null pivot rows” (null pivot row detection)
      MUMPS_CNTL(1) = MUMPSSettings::pivtol; // relative threshold for numerical pivoting (1e-6)
      // debug for MUMPS team ICNTL(2) = 6; ICNTL(3) = 6; ICNTL(4) = 6;
   }

   MUMPSSolver::~MUMPSSolver() {
      this->workspace.job = MUMPSSolver::JOB_END;
      dmumps_c(&this->workspace);
   }

   void MUMPSSolver::initialize_memory() {
      this->workspace.n = static_cast<int>(this->linear_system.dimension);
      this->workspace.nnz = static_cast<int>(this->linear_system.number_nonzeros);
   }

   void MUMPSSolver::do_symbolic_analysis() {
      assert(!this->analysis_performed);

      this->workspace.job = MUMPSSolver::JOB_ANALYSIS;
      // connect the local sparsity with the pointers in the workspace
      this->workspace.irn = this->linear_system.matrix_row_indices.data();
      this->workspace.jcn = this->linear_system.matrix_column_indices.data();
      dmumps_c(&this->workspace);
      if (MUMPS_INFO(1) < 0) {
         throw std::runtime_error("The MUMPS analysis failed");
      }
      if (this->use_scaling) {
         MUMPS_ICNTL(8) = 8; // recompute scaling before factorization
      }
      this->analysis_performed = true;
   }

   void MUMPSSolver::factorize(bool /*is_matrix_positive_definite*/) {
      assert(this->analysis_performed);

      this->workspace.job = MUMPSSolver::JOB_FACTORIZATION;
      this->workspace.a = this->linear_system.matrix_values.data();
      this->factorization_aborted_singular = false;
      size_t number_factorization_failures = 0;
      while (true) {
         dmumps_c(&this->workspace);
         if (MUMPS_INFO(1) == -8 || MUMPS_INFO(1) == -9) { // workspace too small
            if (number_factorization_failures >= MUMPSSettings::max_number_factorization_failures ||
                  MUMPS_ICNTL(14) >= MUMPSSettings::max_memory_increase_percent) {
               throw std::runtime_error("The MUMPS factorization failed (workspace too small)");
            }
            // increase the workspace size (kept for subsequent factorizations) and retry
            MUMPS_ICNTL(14) = std::min(MUMPS_ICNTL(14) * MUMPSSettings::mem_percent_increase,
               MUMPSSettings::max_memory_increase_percent);
            ++number_factorization_failures;
         }
         else if (MUMPS_INFO(1) == -10) {
            // singular matrix: reported via matrix_is_singular() and the inertia (all eigenvalues deemed zero)
            DEBUG << "MUMPS detected a numerically singular matrix\n";
            this->factorization_aborted_singular = true;
            break;
         }
         else if (MUMPS_INFO(1) < 0) {
            throw std::runtime_error("The MUMPS factorization failed");
         }
         else {
            break;
         }
      }
      this->factorization_performed = true;
   }

   void MUMPSSolver::solve(const double* rhs, double* solution, size_t number_of_rhs) {
      assert(this->factorization_performed);

      const size_t dimension = static_cast<size_t>(this->workspace.n);
      if (this->factorization_aborted_singular) {
         // no factors are available: the solution is undefined (the caller handles the singularity)
         view(solution, number_of_rhs * dimension).fill(std::numeric_limits<double>::quiet_NaN());
         return;
      }
      // copy the rhs block into the solution block (overwritten by MUMPS); both are column-major with leading dimension n
      view(solution, number_of_rhs * dimension) = view(rhs, number_of_rhs * dimension);

      this->workspace.nrhs = static_cast<int>(number_of_rhs);
      this->workspace.lrhs = this->workspace.n;
      this->workspace.rhs = solution;
      this->workspace.job = MUMPSSolver::JOB_SOLVE;
      dmumps_c(&this->workspace);
      if (MUMPS_INFO(1) < 0) {
         throw std::runtime_error("The MUMPS solve failed");
      }
   }

   std::optional<size_t> MUMPSSolver::number_delayed_pivots() const {
      if (this->factorization_aborted_singular) {
         return std::nullopt;
      }
      return static_cast<size_t>(MUMPS_INFOG(13)); // total number of delayed pivots
   }

   Inertia MUMPSSolver::get_inertia() const {
      // rank = number_positive_eigenvalues + number_negative_eigenvalues
      // n = rank + number_zero_eigenvalues
      const size_t number_negative_eigenvalues = this->number_negative_eigenvalues();
      const size_t number_zero_eigenvalues = this->number_zero_eigenvalues();
      const size_t number_positive_eigenvalues = static_cast<size_t>(this->workspace.n) - (number_negative_eigenvalues + number_zero_eigenvalues);
      return {number_positive_eigenvalues, number_negative_eigenvalues, number_zero_eigenvalues};
   }

   size_t MUMPSSolver::number_negative_eigenvalues() const {
      if (this->factorization_aborted_singular) {
         return 0;
      }
      return static_cast<size_t>(MUMPS_INFOG(12));
   }

   size_t MUMPSSolver::number_zero_eigenvalues() const {
      if (this->factorization_aborted_singular) {
         return static_cast<size_t>(this->workspace.n);
      }
      return static_cast<size_t>(MUMPS_INFOG(28));
   }

   bool MUMPSSolver::matrix_is_singular() const {
      return (this->number_zero_eigenvalues() > 0);
   }

   size_t MUMPSSolver::rank() const {
      return static_cast<size_t>(this->workspace.n) - this->number_zero_eigenvalues();
   }

   bool MUMPSSolver::increase_pivoting_tolerance() {
      const double current_tolerance = MUMPS_CNTL(1);
      const auto new_tolerance = this->pivoting_tolerance_schedule.increase(current_tolerance);
      if (!new_tolerance.has_value()) {
         return false;
      }
      MUMPS_CNTL(1) = *new_tolerance;
      DEBUG << "Increasing the MUMPS pivoting tolerance from " << current_tolerance << " to " << MUMPS_CNTL(1) << '\n';
      return true;
   }


   LinearSystem& MUMPSSolver::get_linear_system() {
      return this->linear_system;
   }

   COOLinearSystem& MUMPSSolver::get_coo_linear_system() {
      return this->linear_system;
   }

   // protected member functions

   int& MUMPSSolver::MUMPS_ICNTL(size_t index) {
      // handle the Fortran indexing (starting at 1)
      return this->workspace.icntl[index-1];
   }

   double& MUMPSSolver::MUMPS_CNTL(size_t index) {
      // handle the Fortran indexing (starting at 1)
      return this->workspace.cntl[index-1];
   }

   int MUMPSSolver::MUMPS_INFO(size_t index) const {
      // handle the Fortran indexing (starting at 1)
      return this->workspace.info[index-1];
   }

   int MUMPSSolver::MUMPS_INFOG(size_t index) const {
      // handle the Fortran indexing (starting at 1)
      return this->workspace.infog[index-1];
   }
} // namespace