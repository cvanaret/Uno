// Copyright (c) 2024 Manuel Schaich
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include <stdexcept>
#include "MA27Solver.hpp"
#include "ingredients/subproblems/LagrangeNewtonSubproblem.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/Logger.hpp"
#include "fortran_interface.h"

#define ma27_set_parameters FC_GLOBAL(ma27id, MA27ID)
#define ma27_symbolic_analysis FC_GLOBAL(ma27ad, MA27AD)
#define ma27_numerical_factorization FC_GLOBAL(ma27bd, MA27BD)
#define ma27_solve_linear_system FC_GLOBAL(ma27cd, MA27CD)

extern "C" {
   void ma27_set_parameters(int ICNTL[], double CNTL[]);
   void ma27_symbolic_analysis(int* N, int* NZ, int IRN[], int ICN[], int IW[], int* LIW, int IKEEP[], int IW1[],
         int* NSTEPS, int* IFLAG, int ICNTL[], double CNTL[], int INFO[], double* OPS);
   void ma27_numerical_factorization(int* N, int* NZ, int IRN[], int ICN[], double A[], int* LA, int IW[], int* LIW,
         int IKEEP[], int* NSTEPS, int* MAXFRT, int IW1[], int ICNTL[], double CNTL[], int INFO[]);
   void ma27_solve_linear_system(int* N, double A[], int* LA, int IW[], int* LIW, double W[], int* MAXFRT, double RHS[],
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

   MA27Solver::MA27Solver(size_t number_variables, size_t number_constraints, size_t number_jacobian_nonzeros, size_t number_hessian_nonzeros):
      DirectEqualityQPSolver<size_t, double>(),
      objective_gradient(number_variables),
      constraints(number_constraints),
      constraint_jacobian(number_constraints, number_variables), // TODO construct better
         hessian(number_variables, number_hessian_nonzeros, false, "COO"),
      dimension(number_variables + number_constraints), number_nonzeros(number_hessian_nonzeros + number_jacobian_nonzeros),
      augmented_matrix(this->dimension, this->number_nonzeros, true, "COO"),
      rhs(this->dimension),
      solution(this->dimension),
      iw((2 * this->number_nonzeros + 3 * this->dimension + 1) * 6 / 5), // 20% more than 2*nnz + 3*n + 1
         ikeep(3 * this->dimension), iw1(2 * this->dimension) {
      this->row_indices.reserve(this->number_nonzeros);
      this->column_indices.reserve(this->number_nonzeros);
      // initialization: set the default values of the controlling parameters
      ma27_set_parameters(this->icntl.data(), this->cntl.data());
      // a suitable pivot order is to be chosen automatically
      this->iflag = 0;
      // suppress warning messages
      this->icntl[eICNTL::LP] = 0;
      this->icntl[eICNTL::MP] = 0;
      this->icntl[eICNTL::LDIAG] = 0;
   }

   void MA27Solver::do_symbolic_analysis(const SymmetricMatrix<size_t, double>& matrix) {
      assert(matrix.dimension() <= this->iw1.capacity() && "MA27Solver: the dimension of the matrix is larger than the preallocated size");
      assert(matrix.number_nonzeros() <= this->row_indices.capacity() &&
             "MA27Solver: the number of nonzeros of the matrix is larger than the preallocated size");

      // build the internal matrix representation
      save_matrix_to_local_format(matrix);

      int n = static_cast<int>(matrix.dimension());
      int nnz = static_cast<int>(matrix.number_nonzeros());
      this->dimension = matrix.dimension();
      this->number_nonzeros = matrix.number_nonzeros();

      // symbolic analysis
      int liw = static_cast<int>(this->iw.size());
      ma27_symbolic_analysis(&n, &nnz,                                   /* size info */
            this->row_indices.data(), this->column_indices.data(),                     /* matrix indices */
            this->iw.data(), &liw, this->ikeep.data(), this->iw1.data(),  /* solver workspace */
            &this->nsteps, &this->iflag, this->icntl.data(), this->cntl.data(), this->info.data(), &this->ops);

      // resize the factor by at least INFO(5) (here, 50% more)
      this->factor.resize(static_cast<size_t>(3 * this->info[eINFO::NRLNEC] / 2));

      assert(this->info[eINFO::IFLAG] == eIFLAG::SUCCESS && "MA27: the symbolic analysis failed");
      if (this->info[eINFO::IFLAG] != eIFLAG::SUCCESS) {
         WARNING << "MA27 has issued a warning: IFLAG = " << this->info[eINFO::IFLAG] << " additional info, IERROR = " << this->info[eINFO::IERROR] << '\n';
      }
   }

   void MA27Solver::do_numerical_factorization([[maybe_unused]] const SymmetricMatrix<size_t, double>& matrix) {
      assert(matrix.dimension() <= this->iw1.capacity() && "MA27Solver: the dimension of the matrix is larger than the preallocated size");
      assert(this->number_nonzeros == matrix.number_nonzeros() && "MA27Solver: the numbers of nonzeros do not match");

      // initialize factor with the entries of the matrix. It will be modified by MA27BD
      std::copy(matrix.data_pointer(), matrix.data_pointer() + matrix.number_nonzeros(), this->factor.begin());

      int n = static_cast<int>(matrix.dimension());
      int nnz = static_cast<int>(matrix.number_nonzeros());

      // numerical factorization
      // may fail because of insufficient space. In this case, more memory is allocated and the factorization tried again
      bool successful_factorization = false;
      size_t attempt = 0;
      while (!successful_factorization) {
         attempt++;
         if (this->number_factorization_attempts < attempt) {
            throw std::runtime_error("MA27 reached the maximum number of factorization attempts");
         }

         int la = static_cast<int>(this->factor.size());
         int liw = static_cast<int>(this->iw.size());
         ma27_numerical_factorization(&n, &nnz, this->row_indices.data(), this->column_indices.data(), this->factor.data(), &la, this->iw.data(), &liw,
            this->ikeep.data(), &this->nsteps, &this->maxfrt, this->iw1.data(), this->icntl.data(), this->cntl.data(), this->info.data());
         successful_factorization = true;

         if (this->info[eINFO::IFLAG] == eIFLAG::INSUFFICIENTINTEGER) {
            INFO << "MA27: insufficient integer workspace, resizing and retrying. \n";
            // increase the size of iw
            this->iw.resize(static_cast<size_t>(this->info[eINFO::IERROR]));
            successful_factorization = false;
         }
         if (this->info[eINFO::IFLAG] == eIFLAG::INSUFFICIENTREAL) {
            INFO << "MA27: insufficient real workspace, resizing and retrying. \n";
            // increase the size of factor
            this->factor.resize(static_cast<size_t>(this->info[eINFO::IERROR]));
            successful_factorization = false;
         }
      }
      this->w.resize(static_cast<size_t>(this->maxfrt));
      this->check_factorization_status();
   }

   SubproblemStatus MA27Solver::solve_equality_constrained_QP(Statistics& statistics, LagrangeNewtonSubproblem& subproblem,
         const Vector<double>& /*initial_point*/, Vector<double>& direction_primals, Multipliers& direction_multipliers, double& /*subproblem_objective*/,
         WarmstartInformation& warmstart_information) {
      // set up the augmented system
      subproblem.assemble_augmented_matrix(statistics, this->objective_gradient, this->constraints, this->constraint_jacobian, this->hessian,
         this->augmented_matrix, *this, warmstart_information);
      subproblem.assemble_augmented_rhs(this->objective_gradient, this->constraints, this->constraint_jacobian, this->rhs, warmstart_information);
      // solve the augmented system
      this->solve_indefinite_linear_system();
      // form the primal-dual direction (note the minus sign for the multipliers)
      direction_primals = view(this->solution, 0, subproblem.number_variables);
      direction_multipliers.constraints = view(-this->solution, subproblem.number_variables, subproblem.number_variables + subproblem.number_constraints);
      return SubproblemStatus::OPTIMAL; // TODO
   }

   void MA27Solver::solve_indefinite_linear_system() {
      int la = static_cast<int>(this->factor.size());
      int liw = static_cast<int>(this->iw.size());

      // copy rhs into solution (overwritten by MA27)
      this->solution = this->rhs;
      int n = static_cast<int>(this->augmented_matrix.dimension());

      ma27_solve_linear_system(&n, this->factor.data(), &la, this->iw.data(), &liw, this->w.data(), &this->maxfrt, this->solution.data(),
         this->iw1.data(), &this->nsteps, this->icntl.data(), this->info.data());

      assert(this->info[eINFO::IFLAG] == eIFLAG::SUCCESS && "MA27: the linear solve failed");
      if (this->info[eINFO::IFLAG] != eIFLAG::SUCCESS) {
         WARNING << "MA27 has issued a warning: IFLAG = " << this->info[eINFO::IFLAG] << " additional info, IERROR = " << this->info[eINFO::IERROR] << '\n';
      }
   }

   std::tuple<size_t, size_t, size_t> MA27Solver::get_inertia() const {
      // rank = number_positive_eigenvalues + number_negative_eigenvalues
      // n = rank + number_zero_eigenvalues
      const size_t rankA = rank();
      const size_t num_negative_eigenvalues = number_negative_eigenvalues();
      const size_t num_positive_eigenvalues = rankA - num_negative_eigenvalues;
      const size_t num_zero_eigenvalues = this->dimension - rankA;
      return std::make_tuple(num_positive_eigenvalues, num_negative_eigenvalues, num_zero_eigenvalues);
   }

   size_t MA27Solver::number_negative_eigenvalues() const {
      return static_cast<size_t>(this->info[eINFO::NEIG]);
   }

   bool MA27Solver::matrix_is_singular() const {
      return (this->info[eINFO::IFLAG] == eIFLAG::SINGULAR || this->info[eINFO::IFLAG] == eIFLAG::RANK_DEFICIENT);
   }

   size_t MA27Solver::rank() const {
      return (this->info[eINFO::IFLAG] == eIFLAG::RANK_DEFICIENT) ? static_cast<size_t>(this->info[eINFO::IERROR]) : this->dimension;
   }

   void MA27Solver::save_matrix_to_local_format(const SymmetricMatrix<size_t, double>& matrix) {
      // build the internal matrix representation
      this->row_indices.clear();
      this->column_indices.clear();
      this->factor.clear();
      constexpr auto fortran_shift = 1;
      for (const auto [row_index, column_index, element]: matrix) {
         this->row_indices.emplace_back(static_cast<int>(row_index + fortran_shift));
         this->column_indices.emplace_back(static_cast<int>(column_index + fortran_shift));
         this->factor.emplace_back(element);
      }
   }

   void MA27Solver::check_factorization_status() {
      switch (this->info[eINFO::IFLAG]) {
         case NSTEPS:
            WARNING << "MA27BD: Value of NSTEPS outside the range 1 ≤ NSTEPS ≤ N" << '\n';
            break;
         case PIVOTSIGN:
            WARNING << "MA27BD: A change of sign of pivots has been detected when U was negative. Detected at pivot step " << this->info[eINFO::IERROR]
                    << '\n';
            break;
         case SINGULAR:
            DEBUG << "MA27BD: Matrix is singular. Singularity detected during pivot step " << this->info[eINFO::IERROR] << '\n';
            break;
         case NZOUTOFRANGE:
            WARNING << "MA27BD: Value of NZ out of range. NZ < 0." << '\n';
            break;
         case NOUTOFRANGE:
            WARNING << "MA27BD: Value of N out of range. N < 1." << '\n';
            break;
         case IDXOUTOFRANGE:
            WARNING << "MA27BD: Index (in IRN or ICN) out of range. " << this->info[eINFO::IERROR] << " indices affected." << '\n';
            break;
         case FALSEDEFINITENESS:
            WARNING << "MA27BD: Matrix was supposed to be definite, but pivots have different signs when factorizing. Detected "
                    << this->info[eINFO::IERROR] << " sign changes." << '\n';
            break;
         case RANK_DEFICIENT:
            DEBUG << "MA27BD: Matrix is rank deficient. Rank: " << this->info[eINFO::IERROR] << " whereas dimension " << this->dimension << '\n';
            break;
      }
   }
} // namespace
