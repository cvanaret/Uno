// Copyright (c) 2024-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include <stdexcept>
#include "MUMPSSolver.hpp"
#include "linear_algebra/Vector.hpp"
#include "options/Options.hpp"
#include "symbolic/Range.hpp"
#include "tools/Logger.hpp"
#if defined(HAS_MPI) && defined(MUMPS_PARALLEL)
#include "mpi.h"
#endif

#define USE_COMM_WORLD (-987654)

namespace uno {
   MUMPSSolver::MUMPSSolver(): DirectSymmetricIndefiniteLinearSolver() {
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
      ICNTL(1) = -1; // output stream for error messages (off)
      ICNTL(2) = -1; // output stream for diagnostic printing (off)
      ICNTL(3) = -1; // output stream for global information (off)
      ICNTL(4) = 0; // level of printing for error, warning, and diagnostic message (off)

      ICNTL(6) = MUMPSSettings::permuting_scaling; // column permutation (none)
      ICNTL(7) = MUMPSSettings::pivot_order; // pivot order (AMD)
      ICNTL(8) = MUMPSSettings::scaling; // scaling strategy (diagonal scaling during factorization)
      ICNTL(10) = 0; // iterative refinement (none)
      ICNTL(13) = 1; // parallelism of the root nod (ScaLAPACK not used, partly recover parallelism of the root node)
      ICNTL(14) = MUMPSSettings::mem_percent; // percentage increase in the estimated working space (35)
      ICNTL(24) = 1; // controls the detection of “null pivot rows” (null pivot row detection)
      CNTL(1) = MUMPSSettings::pivtol; // relative threshold for numerical pivoting (1e-6)
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
      if (INFO(1) < 0) {
         throw std::runtime_error("The MUMPS analysis failed");
      }
      ICNTL(8) = 8; // recompute scaling before factorization
      this->analysis_performed = true;
   }

   void MUMPSSolver::do_numerical_factorization(bool /*is_matrix_positive_definite*/) {
      assert(this->analysis_performed);

      this->workspace.job = MUMPSSolver::JOB_FACTORIZATION;
      this->workspace.a = this->linear_system.matrix_values.data();
      bool success = false;
      while (!success) {
         dmumps_c(&this->workspace);
         if (INFO(1) == -8 || INFO(1) == -9) { // workspace too small
            if (this->number_factorization_failures >= this->max_number_factorization_failures) {
               this->number_factorization_failures = 0;
               throw std::runtime_error("The MUMPS factorization failed (workspace too small)");
            }
            // increase the workspace size and retry
            ICNTL(14) = ICNTL(14) * MUMPSSettings::mem_percent_increase;
            ++this->number_factorization_failures;
         }
         else if (INFO(1) == -10) {
            // singular matrix, should be caught by the calling code via the inertia
            DEBUG << "MUMPS detected a numerically singular matrix\n";
            success = true;
         }
         else if (INFO(1) < 0) {
            throw std::runtime_error("The MUMPS factorization failed");
         }
         else {
            success = true;
         }
      }
      this->factorization_performed = true;
   }

   void MUMPSSolver::solve_indefinite_system(double* solution) {
      assert(this->factorization_performed);

      // copy rhs into solution (overwritten by MUMPS)
      const size_t dimension = static_cast<size_t>(this->workspace.n);
      view(solution, dimension) = this->linear_system.rhs.view();

      this->workspace.nrhs = 1;
      this->workspace.rhs = solution;
      this->workspace.job = MUMPSSolver::JOB_SOLVE;
      dmumps_c(&this->workspace);
      if (INFO(1) < 0) {
         throw std::runtime_error("The MUMPS solve failed");
      }
   }

   void MUMPSSolver::solve_indefinite_system(const double* rhs, double* solution, size_t number_of_rhs) {
      assert(this->factorization_performed);

      // copy the rhs block into the solution block (overwritten by MUMPS); both are column-major,
      // with leading dimension n (lrhs)
      const size_t dimension = static_cast<size_t>(this->workspace.n);
      view(solution, number_of_rhs * dimension) = view(rhs, number_of_rhs * dimension);

      this->workspace.nrhs = static_cast<int>(number_of_rhs);
      this->workspace.lrhs = this->workspace.n;
      this->workspace.rhs = solution;
      this->workspace.job = MUMPSSolver::JOB_SOLVE;
      dmumps_c(&this->workspace);
      if (INFO(1) < 0) {
         throw std::runtime_error("The MUMPS solve failed");
      }
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
      return static_cast<size_t>(INFOG(12));
   }

   size_t MUMPSSolver::number_zero_eigenvalues() const {
      return static_cast<size_t>(INFOG(28));
   }

   bool MUMPSSolver::matrix_is_singular() const {
      return (this->number_zero_eigenvalues() > 0);
   }

   size_t MUMPSSolver::rank() const {
      return static_cast<size_t>(this->workspace.n) - this->number_zero_eigenvalues();
   }

   LinearSystem& MUMPSSolver::get_linear_system() {
      return this->linear_system;
   }

   COOLinearSystem& MUMPSSolver::get_coo_linear_system() {
      return this->linear_system;
   }

   // protected member function

   int& MUMPSSolver::ICNTL(size_t index) {
      // handle the Fortran indexing (starting at 1)
      return this->workspace.icntl[index-1];
   }

   double& MUMPSSolver::CNTL(size_t index) {
      // handle the Fortran indexing (starting at 1)
      return this->workspace.cntl[index-1];
   }

   int MUMPSSolver::INFO(size_t index) const {
      // handle the Fortran indexing (starting at 1)
      return this->workspace.info[index-1];
   }

   int MUMPSSolver::INFOG(size_t index) const {
      // handle the Fortran indexing (starting at 1)
      return this->workspace.infog[index-1];
   }
} // namespace