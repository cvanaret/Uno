// Copyright (c) 2024-2026 Manuel Schaich, Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include "MA27Solver.hpp"
#include "tools/Logger.hpp"
#ifdef HSL_RUNTIME_LOADING
#include "ingredients/subproblem_solvers/HSL/HSLLoader.hpp"
// route the calls through the runtime-resolved function pointers
#define LIBHSL_version uno::hsl_libhsl_version
#define MA27_set_default_parameters uno::hsl_ma27id
#define MA27_symbolic_analysis uno::hsl_ma27ad
#define MA27_numerical_factorization uno::hsl_ma27bd
#define MA27_linear_solve uno::hsl_ma27cd
#else
#include "fortran_interface.h"
#define MA27_set_default_parameters FC_GLOBAL(ma27id, MA27ID)
#define MA27_symbolic_analysis FC_GLOBAL(ma27ad, MA27AD)
#define MA27_numerical_factorization FC_GLOBAL(ma27bd, MA27BD)
#define MA27_linear_solve FC_GLOBAL(ma27cd, MA27CD)

extern "C" {
   void LIBHSL_version(int *major, int *minor, int *patch);

   void MA27_set_default_parameters(int ICNTL[], double CNTL[]);

   void MA27_symbolic_analysis(int* N, int* NZ, int IRN[], int ICN[], int IW[], int* LIW, int IKEEP[], int IW1[],
      int* NSTEPS, int* IFLAG, int ICNTL[], double CNTL[], int INFO[], double* OPS);

   void MA27_numerical_factorization(int* N, int* NZ, int IRN[], int ICN[], double A[], int* LA, int IW[], int* LIW,
      int IKEEP[], int* NSTEPS, int* MAXFRT, int IW1[], int ICNTL[], double CNTL[], int INFO[]);

   void MA27_linear_solve(int* N, double A[], int* LA, int IW[], int* LIW, double W[], int* MAXFRT, double RHS[],
      int IW1[], int* NSTEPS, int ICNTL[], int INFO[]);
}
#endif

namespace uno {
   enum eIFLAG {
      NSTEPS = -7, // Value of NSTEPS outside the range 1 ≤ NSTEPS ≤ N (MA27B/BD entry).
      PIVOTSIGN = -6, // A change of sign of pivots has been detected when U was negative. INFO(2) is set to the pivot step at which the change was detected. (MA27B/BD entry only)
      SINGULAR = -5, // Matrix is singular (MA27B/BD entry only). INFO(2) is set to the pivot step at which singularity was detected
      INSUFFICIENTREAL = -4, // Failure due to insufficient space allocated to array A (MA27B/BD entry only). INFO(2) is set to a value that may suffice.
      INSUFFICIENTINTEGER = -3, // Failure due to insufficient space allocated to array IW (MA27A/AD and MA27B/BD entries). INFO(2) is set to a value that may suffice.
      NZOUTOFRANGE = -2, // Value of NZ out of range. NZ < 0. (MA27A/AD and MA27B/BD entries)
      NOUTOFRANGE = -1, // Value of N out of range. N < 1. (MA27A/AD and MA27B/BD entries).
      SUCCESS = 0, // Successful completion.
      IDXOUTOFRANGE = 1, // ndex (in IRN or ICN) out of range. Action taken by subroutine is to ignore any such entries and continue (MA27A/AD and MA27B/BD entries). INFO(2) is set to the number of faulty entries. Details of the first ten are printed on unit ICNTL(2).
      FALSEDEFINITENESS, // Pivots have different signs when factorizing a supposedly definite matrix (when the value of U in CNTL(1) is zero) (MA27B/BD entry only). INFO(2) is set to the number of sign changes. Note that this warning will overwrite an INFO(1)=1 warning. Details of the first ten are printed on unit ICNTL(2).
      RANK_DEFICIENT, // Matrix is rank deficient. In this case, a decomposition will still have been produced which will enable the subsequent solution of consistent equations (MA27B/BD entry only). INFO(2) will be set to the rank of the matrix. Note that this warning will overwrite an INFO(1)=1 or INFO(1)=2 warning.
   };

   MA27Solver::MA27Solver(): DirectSymmetricIndefiniteLinearSolver() {
#ifdef HSL_RUNTIME_LOADING
      if (!ma27_symbols_available()) {
         throw std::runtime_error("Uno: the MA27 solver was requested but the HSL library could not be loaded at "
            "runtime (set the UNO_HSL_LIBRARY environment variable to point at libhsl)");
      }
#endif
#if defined(HAS_HSL)
      int major, minor, patch;
      LIBHSL_version(&major, &minor, &patch);
      INFO << "Running MA27 v" << major << "." << minor << "." << patch << '\n';
#else
      INFO << "Running MA27 v1.0.0\n";
#endif
      // initialization: set the default values of the controlling parameters
      MA27_set_default_parameters(this->workspace.icntl.data(), this->workspace.cntl.data());
      MA27_ICNTL(1) = 0; // no warning messages
      MA27_ICNTL(2) = 0; // no diagnostic printing
      MA27_ICNTL(3) = 0; // no diagnostic printing
      MA27_CNTL(1) = MA27Settings::pivoting_threshold;
      this->workspace.iflag = 0; // a suitable pivot order is to be chosen automatically
   }

   void MA27Solver::initialize_memory() {
      this->workspace.n = static_cast<int>(this->linear_system.dimension);
      this->workspace.nnz = static_cast<int>(this->linear_system.number_nonzeros);
      // 100% more than 2*nnz + 3*n + 1
      const size_t liw = 2 * (2 * this->linear_system.number_nonzeros + 3 * this->linear_system.dimension + 1);
      this->workspace.liw = static_cast<int>(liw);
      this->workspace.iw.resize(liw);
      this->workspace.ikeep.resize(3 * this->linear_system.dimension);
      this->workspace.iw1.resize(2 * this->linear_system.dimension);
   }

   void MA27Solver::do_symbolic_analysis() {
      assert(!this->analysis_performed);

      MA27_symbolic_analysis(&this->workspace.n, &this->workspace.nnz, /* size info */
         this->linear_system.matrix_row_indices.data(), this->linear_system.matrix_column_indices.data(),
         this->workspace.iw.data(), &this->workspace.liw, this->workspace.ikeep.data(), this->workspace.iw1.data(),  /* solver workspace */
         &this->workspace.nsteps, &this->workspace.iflag, this->workspace.icntl.data(), this->workspace.cntl.data(),
         this->workspace.info.data(), &this->workspace.ops);
      if (MA27_INFO(1) != 0) {
         throw std::runtime_error("MA27: the symbolic analysis failed");
      }
      // resize IW
      const int nirnec = MA27_INFO(6);
      this->workspace.liw = MA27Settings::liw_init_factor * nirnec;
      this->workspace.iw.reserve(static_cast<size_t>(this->workspace.liw));
      // resize factor (A)
      const int nrlnec = MA27_INFO(5);
      this->workspace.la = MA27Settings::la_init_factor * nrlnec;
      this->workspace.factor.resize(static_cast<size_t>(this->workspace.la));

      this->analysis_performed = true;
   }

   void MA27Solver::do_numerical_factorization(bool /*is_matrix_positive_definite*/) {
      assert(this->analysis_performed);

      // initialize factor with the entries of the matrix. It will be modified by MA27BD
      std::copy_n(this->linear_system.matrix_values.data(), this->workspace.nnz, this->workspace.factor.begin());

      // numerical factorization
      // may fail because of insufficient space. In this case, more memory is allocated and the factorization tried again
      bool factorization_done = false;
      size_t attempt = 0;
      while (not factorization_done) {
         ++attempt;
         if (this->workspace.number_factorization_attempts < attempt) {
            throw std::runtime_error("MA27 reached the maximum number of factorization attempts");
         }

         int la = static_cast<int>(this->workspace.factor.size());
         int liw = static_cast<int>(this->workspace.iw.size());
         MA27_numerical_factorization(&this->workspace.n, &this->workspace.nnz, this->linear_system.matrix_row_indices.data(),
            this->linear_system.matrix_column_indices.data(), this->workspace.factor.data(), &la, this->workspace.iw.data(), &liw,
            this->workspace.ikeep.data(), &this->workspace.nsteps, &this->workspace.maxfrt, this->workspace.iw1.data(),
            this->workspace.icntl.data(), this->workspace.cntl.data(), this->workspace.info.data());
         factorization_done = true;

         if (MA27_INFO(1) == eIFLAG::INSUFFICIENTINTEGER) {
            DEBUG << "MA27: insufficient integer workspace, resizing and retrying. \n";
            // increase the size of iw by 50%
            this->workspace.iw.resize(static_cast<size_t>(3 * this->workspace.info[1] / 2));
            factorization_done = false;
         }
         if (MA27_INFO(1) == eIFLAG::INSUFFICIENTREAL) {
            DEBUG << "MA27: insufficient real workspace, resizing and retrying. \n";
            // increase the size of factor by 50%
            this->workspace.factor.resize(static_cast<size_t>(3 * this->workspace.info[1] / 2));
            factorization_done = false;
         }
      }
      this->workspace.w.resize(static_cast<size_t>(this->workspace.maxfrt));
      this->check_factorization_status();
      this->factorization_performed = true;
   }

   void MA27Solver::solve_indefinite_system(double* solution) {
      assert(this->factorization_performed);

      // copy rhs into solution (overwritten by MA27)
      const size_t dimension = static_cast<size_t>(this->workspace.n);
      view(solution, dimension) = this->linear_system.rhs.view();

      MA27_linear_solve(&this->workspace.n, this->workspace.factor.data(), &this->workspace.la, this->workspace.iw.data(),
         &this->workspace.liw, this->workspace.w.data(), &this->workspace.maxfrt, solution, this->workspace.iw1.data(),
         &this->workspace.nsteps, this->workspace.icntl.data(), this->workspace.info.data());

      assert(MA27_INFO(1) == 0 && "MA27: the linear solve failed");
      if (MA27_INFO(1) != 0) {
         WARNING << "MA27 has issued a warning: IFLAG = " << MA27_INFO(1) << " additional info, IERROR = "
            << this->workspace.info[1] << '\n';
      }
   }

   Inertia MA27Solver::get_inertia() const {
      // rank = number_positive_eigenvalues + number_negative_eigenvalues
      // n = rank + number_zero_eigenvalues
      const size_t rank = this->rank();
      const size_t number_negative_eigenvalues = this->number_negative_eigenvalues();
      const size_t number_positive_eigenvalues = rank - number_negative_eigenvalues;
      const size_t number_zero_eigenvalues = static_cast<size_t>(this->workspace.n) - rank;
      return {number_positive_eigenvalues, number_negative_eigenvalues, number_zero_eigenvalues};
   }

   size_t MA27Solver::number_negative_eigenvalues() const {
      return static_cast<size_t>(MA27_INFO(15));
   }

   bool MA27Solver::matrix_is_singular() const {
      return (MA27_INFO(1) == eIFLAG::SINGULAR || MA27_INFO(1) == eIFLAG::RANK_DEFICIENT);
   }

   size_t MA27Solver::rank() const {
      return (MA27_INFO(1) == eIFLAG::RANK_DEFICIENT) ?
         static_cast<size_t>(this->workspace.info[1]) :
         static_cast<size_t>(this->workspace.n);
   }

   LinearSystem& MA27Solver::get_linear_system() {
      return this->linear_system;
   }

   COOLinearSystem& MA27Solver::get_coo_linear_system() {
      return this->linear_system;
   }

   // protected member functions

   void MA27Solver::check_factorization_status() const {
      switch (MA27_INFO(1)) {
         case NSTEPS:
            WARNING << "MA27BD: Value of NSTEPS outside the range 1 ≤ NSTEPS ≤ N" << '\n';
            break;
         case PIVOTSIGN:
            WARNING << "MA27BD: A change of sign of pivots has been detected when U was negative. Detected at pivot step "
               << this->workspace.info[1] << '\n';
            break;
         case SINGULAR:
            DEBUG << "MA27BD: Matrix is singular. Singularity detected during pivot step " << this->workspace.info[1] << '\n';
            break;
         case NZOUTOFRANGE:
            WARNING << "MA27BD: Value of NZ out of range. NZ < 0." << '\n';
            break;
         case NOUTOFRANGE:
            WARNING << "MA27BD: Value of N out of range. N < 1." << '\n';
            break;
         case IDXOUTOFRANGE:
            WARNING << "MA27BD: Index (in IRN or ICN) out of range. " << this->workspace.info[1] << " indices affected." << '\n';
            break;
         case FALSEDEFINITENESS:
            WARNING << "MA27BD: Matrix was supposed to be definite, but pivots have different signs when factorizing. Detected "
                    << this->workspace.info[1] << " sign changes." << '\n';
            break;
         case RANK_DEFICIENT:
            DEBUG << "MA27BD: Matrix is rank deficient. Rank: " << this->workspace.info[1] << " whereas dimension "
               << this->workspace.n << '\n';
            break;
      }
   }

   int& MA27Solver::MA27_ICNTL(size_t index) {
      // handle the Fortran indexing (starting at 1)
      return this->workspace.icntl[index-1];
   }

   double& MA27Solver::MA27_CNTL(size_t index) {
      // handle the Fortran indexing (starting at 1)
      return this->workspace.cntl[index-1];
   }

   int MA27Solver::MA27_INFO(size_t index) const {
      // handle the Fortran indexing (starting at 1)
      return this->workspace.info[index-1];
   }
} // namespace