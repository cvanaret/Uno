// Copyright (c) 2024 Manuel Schaich
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include "MA27Solver.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Direction.hpp"
#include "tools/Logger.hpp"
#include "fortran_interface.h"

#define MA27_set_default_parameters FC_GLOBAL(ma27id, MA27ID)
#define MA27_symbolic_analysis FC_GLOBAL(ma27ad, MA27AD)
#define MA27_numerical_factorization FC_GLOBAL(ma27bd, MA27BD)
#define MA27_linear_solve FC_GLOBAL(ma27cd, MA27CD)

extern "C" {
   void MA27_set_default_parameters(int ICNTL[], double CNTL[]);

   void MA27_symbolic_analysis(int* N, int* NZ, int IRN[], int ICN[], int IW[], int* LIW, int IKEEP[], int IW1[],
      int* NSTEPS, int* IFLAG, int ICNTL[], double CNTL[], int INFO[], double* OPS);

   void MA27_numerical_factorization(int* N, int* NZ, int IRN[], int ICN[], double A[], int* LA, int IW[], int* LIW,
      int IKEEP[], int* NSTEPS, int* MAXFRT, int IW1[], int ICNTL[], double CNTL[], int INFO[]);

   void MA27_linear_solve(int* N, double A[], int* LA, int IW[], int* LIW, double W[], int* MAXFRT, double RHS[],
      int IW1[], int* NSTEPS, int ICNTL[], int INFO[]);
}

namespace uno {
   enum eICNTL {
      LP = 0, // Used by the subroutines as the output stream for error messages. If it is set to zero these messages will be suppressed. The default value is 6.
      MP, // Used by the subroutines as the output stream for diagnostic printing and for warning messages. If it is set to zero then messages are suppressed. The default value is 6.
      LDIAG, // Used by the subroutines to control diagnostic printing. If ICNTL(3) is equal to zero (the default), no diagnostic printing will be produced, a value of 1 will print scalar parameters (both in argument lists and in the control and information arrays) and a few entries of array parameters on entry and successful exit from each subroutine while ICNTL(3) equal to 2 will print all parameter values on entry and successful exit.
      /* The entries ICNTL(4) to ICNTL(25) are not of interest to the general user and are discussed more fully by Duff and Reid (AERE R-10533, 1982) under the internal names IOVFLO, NEMIN and IFRLVL
      */
      IOVFLO,
      NEMIN,
      IFRLVL1,
      IFRLVL2,
      IFRLVL3,
      IFRLVL4,
      IFRLVL5,
      IFRLVL6,
      IFRLVL7,
      IFRLVL8,
      IFRLVL9,
      IFRLVL10,
      IFRLVL11,
      IFRLVL12,
      IFRLVL13,
      IFRLVL14,
      IFRLVL15,
      IFRLVL16,
      IFRLVL17,
      IFRLVL18,
      IFRLVL19,
      IFRLVL20,
      UNUSED_ICNTL1,
      UNUSED_ICNTL2,
      UNUSED_ICNTL3,
      UNUSED_ICNTL4,
      UNUSED_ICNTL5,
   };

   enum eCNTL {
      U = 0, // Used by the subroutine to control numerical pivoting. Values greater than 0.5 are treated as 0.5 and less than –0.5 as –0.5. Its default value is 0.1. If U is positive, numerical pivoting will be performed. If U is non-positive, no pivoting will be performed, the subroutine will fail if a zero pivot is encountered, and a flag (see section 2.3) will be set if not all pivots are of the same sign; the factorization will continue after a sign change is detected if U is zero but will exit immediately if U is less than zero. If the system is definite, then setting U to zero will decrease the factorization time while still providing a stable decomposition. For problems requiring greater than average numerical care a higher value than the default would be advisable.
      FRATIO, // Given the default value of 1.0 by MA27I/ID. If MA27A/AD encounters a row of the reduced matrix with a proportion of entries greater than FRATIO, the row is treated as full. FRATIO is not altered by MA27.
      PIVTOL, // Given the default value of 0.0 by MA27I/ID. MA27B/BD will not accept an entry with absolute value less than PIVTOL as a 1×1 pivot or the off-diagonal entry of a 2×2 pivot. PIVTOL is not altered by MA27.
      UNUSED_CNTL1,
      UNUSED_CNTL2,
   };

   enum eINFO {
      IFLAG = 0, // An error flag. A value of zero indicates that the subroutine has performed successfully.
      IERROR, // Provides supplementary information when there is an error.
      NRLTOT, // Gives the total amount of REAL words required for a successful completion of MA27B/BD without the need for data compression provided no numerical pivoting is performed. The actual amount required may be higher because of numerical pivoting, but probably not by more than 3%.
      NIRTOT, // Gives the total amount of INTEGER words required for a successful completion of MA27B/BD without the need for data compression provided no numerical pivoting is performed. The actual amount required may be higher because of numerical pivoting, but probably not by more than 3%.
      NRLNEC, // Gives the amount of REAL words required for successful completion of MA27B/BD allowing data compression (see NCMPBR returned in INFO(12)), again provided no numerical pivoting is performed. Numerical pivoting may cause a higher value to be required, but probably not by more than 3%. If storage was conserved by equivalencing IW(1) with IRN(1), NRLNEC and NIRNEC cannot be calculated exactly but instead an upper bound will be returned. Experience has shown that this can overestimate the exact values by 50% although the tightness of the bound is very problem dependent. For example, a tight bound will generally be obtained if there are many more entries in the factors than in the input matrix.
      NIRNEC, // Gives the amount of INTEGER words required for successful completion of MA27B/BD allowing data compression (see NCMPBR returned in INFO(12)), again provided no numerical pivoting is performed. Numerical pivoting may cause a higher value to be required, but probably not by more than 3%. If storage was conserved by equivalencing IW(1) with IRN(1), NRLNEC and NIRNEC cannot be calculated exactly but instead an upper bound will be returned. Experience has shown that this can overestimate the exact values by 50% although the tightness of the bound is very problem dependent. For example, a tight bound will generally be obtained if there are many more entries in the factors than in the input matrix.
      NRLADU, // Gives the number of REAL words required to hold the matrix factors if no numerical pivoting is performed by MA27B/BD. Numerical pivoting may change this slightly.
      NIRADU, // Gives the number of INTEGER words required to hold the matrix factors if no numerical pivoting is performed by MA27B/BD. Numerical pivoting may change this slightly.
      NRLBDU, // Gives the amount of REAL words actually used to hold the factorization.
      NIRBDU, // Gives the amount of INTEGER words actually used to hold the factorization.
      NCMPA, // Holds the number of compresses of the internal data structure performed by MA27A/AD. If this is high (say > 10), the performance of MA27A/AD may be improved by increasing the length of array IW.
      NCMPBR, // Holds the number of compresses of the real data structure required by the factorization. If either of these is high (say > 10), then the speed of the factorization may be increased by allocating more space to the arrays A as appropriate.
      NCMPBI, // Holds the number of compresses of the integer data structure required by the factorization. If either of these is high (say > 10), then the speed of the factorization may be increased by allocating more space to the arrays IW as appropriate.
      NTWO, // Gives the number of 2×2 pivots used during the factorization.
      NEIG, // Gives the number of negative eigenvalues of A.
      UNUSED_INFO1,
      UNUSED_INFO2,
      UNUSED_INFO3,
      UNUSED_INFO4,
      UNUSED_INFO5,
   };

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
      // initialization: set the default values of the controlling parameters
      MA27_set_default_parameters(this->workspace.icntl.data(), this->workspace.cntl.data());
      // a suitable pivot order is to be chosen automatically
      this->workspace.iflag = 0;
      // suppress warning messages
      this->workspace.icntl[eICNTL::LP] = 0;
      this->workspace.icntl[eICNTL::MP] = 0;
      this->workspace.icntl[eICNTL::LDIAG] = 0;
   }

   void MA27Solver::initialize_hessian(const Subproblem& subproblem) {
      this->evaluation_space.initialize_hessian(subproblem);

      // workspace
      const size_t dimension = subproblem.number_variables;
      this->workspace.n = static_cast<int>(dimension);
      this->workspace.nnz = static_cast<int>(this->evaluation_space.number_matrix_nonzeros);
      // 20% more than 2*nnz + 3*n + 1
      this->workspace.iw.resize((2 * this->evaluation_space.number_matrix_nonzeros + 3 * dimension + 1) * 6 / 5);
      this->workspace.ikeep.resize(3 * dimension);
      this->workspace.iw1.resize(2 * dimension);
   }

   void MA27Solver::initialize_augmented_system(const Subproblem& subproblem) {
      this->evaluation_space.initialize_augmented_system(subproblem);

      // workspace
      const size_t dimension = subproblem.number_variables + subproblem.number_constraints;
      this->workspace.n = static_cast<int>(dimension);
      this->workspace.nnz = static_cast<int>(this->evaluation_space.number_matrix_nonzeros);
      // 20% more than 2*nnz + 3*n + 1
      this->workspace.iw.resize((2 * this->evaluation_space.number_matrix_nonzeros + 3 * dimension + 1) * 6 / 5);
      this->workspace.ikeep.resize(3 * dimension);
      this->workspace.iw1.resize(2 * dimension);
   }

   void MA27Solver::do_symbolic_analysis() {
      assert(!this->analysis_performed);

      int liw = static_cast<int>(this->workspace.iw.size());
      MA27_symbolic_analysis(&this->workspace.n, &this->workspace.nnz,              /* size info */
         this->evaluation_space.matrix_row_indices.data(), this->evaluation_space.matrix_column_indices.data(),                     /* matrix indices */
         this->workspace.iw.data(), &liw, this->workspace.ikeep.data(), this->workspace.iw1.data(),  /* solver workspace */
         &this->workspace.nsteps, &this->workspace.iflag, this->workspace.icntl.data(), this->workspace.cntl.data(),
         this->workspace.info.data(), &this->workspace.ops);

      // resize the factor by at least INFO(5) (here, 50% more)
      this->workspace.factor.resize(static_cast<size_t>(3 * this->workspace.info[eINFO::NRLNEC] / 2));

      assert(this->workspace.info[eINFO::IFLAG] == eIFLAG::SUCCESS && "MA27: the symbolic analysis failed");
      if (this->workspace.info[eINFO::IFLAG] != eIFLAG::SUCCESS) {
         WARNING << "MA27 has issued a warning: IFLAG = " << this->workspace.info[eINFO::IFLAG] << " additional info, IERROR = "
            << this->workspace.info[eINFO::IERROR] << '\n';
      }
      this->analysis_performed = true;
   }

   void MA27Solver::do_numerical_factorization(const double* matrix_values) {
      assert(this->analysis_performed);

      // initialize factor with the entries of the matrix. It will be modified by MA27BD
      std::copy_n(matrix_values, this->workspace.nnz, this->workspace.factor.begin());

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
         MA27_numerical_factorization(&this->workspace.n, &this->workspace.nnz, this->evaluation_space.matrix_row_indices.data(),
            this->evaluation_space.matrix_column_indices.data(), this->workspace.factor.data(), &la, this->workspace.iw.data(), &liw,
            this->workspace.ikeep.data(), &this->workspace.nsteps, &this->workspace.maxfrt, this->workspace.iw1.data(),
            this->workspace.icntl.data(), this->workspace.cntl.data(), this->workspace.info.data());
         factorization_done = true;

         if (this->workspace.info[eINFO::IFLAG] == eIFLAG::INSUFFICIENTINTEGER) {
            DEBUG << "MA27: insufficient integer workspace, resizing and retrying. \n";
            // increase the size of iw by 50%
            this->workspace.iw.resize(static_cast<size_t>(3 * this->workspace.info[eINFO::IERROR] / 2));
            factorization_done = false;
         }
         if (this->workspace.info[eINFO::IFLAG] == eIFLAG::INSUFFICIENTREAL) {
            DEBUG << "MA27: insufficient real workspace, resizing and retrying. \n";
            // increase the size of factor by 50%
            this->workspace.factor.resize(static_cast<size_t>(3 * this->workspace.info[eINFO::IERROR] / 2));
            factorization_done = false;
         }
      }
      this->workspace.w.resize(static_cast<size_t>(this->workspace.maxfrt));
      this->check_factorization_status();
      this->factorization_performed = true;
   }

   void MA27Solver::solve_indefinite_system(const Vector<double>& /*matrix_values*/, const Vector<double>& rhs,
         Vector<double>& result) {
      assert(this->factorization_performed);

      int la = static_cast<int>(this->workspace.factor.size());
      int liw = static_cast<int>(this->workspace.iw.size());

      result = rhs;

      MA27_linear_solve(&this->workspace.n, this->workspace.factor.data(), &la, this->workspace.iw.data(), &liw,
         this->workspace.w.data(), &this->workspace.maxfrt, result.data(), this->workspace.iw1.data(), &this->workspace.nsteps,
         this->workspace.icntl.data(), this->workspace.info.data());

      assert(this->workspace.info[eINFO::IFLAG] == eIFLAG::SUCCESS && "MA27: the linear solve failed");
      if (this->workspace.info[eINFO::IFLAG] != eIFLAG::SUCCESS) {
         WARNING << "MA27 has issued a warning: IFLAG = " << this->workspace.info[eINFO::IFLAG] << " additional info, IERROR = "
            << this->workspace.info[eINFO::IERROR] << '\n';
      }
   }

   void MA27Solver::solve_indefinite_system(Statistics& statistics, const Subproblem& subproblem, Direction& direction,
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
      return static_cast<size_t>(this->workspace.info[eINFO::NEIG]);
   }

   bool MA27Solver::matrix_is_singular() const {
      return (this->workspace.info[eINFO::IFLAG] == eIFLAG::SINGULAR || this->workspace.info[eINFO::IFLAG] == eIFLAG::RANK_DEFICIENT);
   }

   size_t MA27Solver::rank() const {
      return (this->workspace.info[eINFO::IFLAG] == eIFLAG::RANK_DEFICIENT) ?
         static_cast<size_t>(this->workspace.info[eINFO::IERROR]) :
         static_cast<size_t>(this->workspace.n);
   }

   EvaluationSpace& MA27Solver::get_evaluation_space() {
      return this->evaluation_space;
   }

   void MA27Solver::check_factorization_status() {
      switch (this->workspace.info[eINFO::IFLAG]) {
         case NSTEPS:
            WARNING << "MA27BD: Value of NSTEPS outside the range 1 ≤ NSTEPS ≤ N" << '\n';
            break;
         case PIVOTSIGN:
            WARNING << "MA27BD: A change of sign of pivots has been detected when U was negative. Detected at pivot step "
               << this->workspace.info[eINFO::IERROR] << '\n';
            break;
         case SINGULAR:
            DEBUG << "MA27BD: Matrix is singular. Singularity detected during pivot step " << this->workspace.info[eINFO::IERROR] << '\n';
            break;
         case NZOUTOFRANGE:
            WARNING << "MA27BD: Value of NZ out of range. NZ < 0." << '\n';
            break;
         case NOUTOFRANGE:
            WARNING << "MA27BD: Value of N out of range. N < 1." << '\n';
            break;
         case IDXOUTOFRANGE:
            WARNING << "MA27BD: Index (in IRN or ICN) out of range. " << this->workspace.info[eINFO::IERROR] << " indices affected." << '\n';
            break;
         case FALSEDEFINITENESS:
            WARNING << "MA27BD: Matrix was supposed to be definite, but pivots have different signs when factorizing. Detected "
                    << this->workspace.info[eINFO::IERROR] << " sign changes." << '\n';
            break;
         case RANK_DEFICIENT:
            DEBUG << "MA27BD: Matrix is rank deficient. Rank: " << this->workspace.info[eINFO::IERROR] << " whereas dimension "
               << this->workspace.n << '\n';
            break;
      }
   }
} // namespace