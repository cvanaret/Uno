// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "MUMPSSolver.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#if defined(HAS_MPI) && defined(MUMPS_PARALLEL)
#include "mpi.h"
#endif

#define USE_COMM_WORLD (-987654)

namespace uno {
   MUMPSSolver::MUMPSSolver(size_t dimension, size_t number_nonzeros) : DirectSymmetricIndefiniteLinearSolver<size_t, double>(dimension) {
      this->row_indices.reserve(number_nonzeros);
      this->column_indices.reserve(number_nonzeros);

      this->mumps_structure.sym = MUMPSSolver::GENERAL_SYMMETRIC;
#if defined(HAS_MPI) && defined(MUMPS_PARALLEL)
      // TODO load number of processes from option file
      this->mumps_structure.par = 1;
#else
      this->mumps_structure.par = 1;
#endif
      this->mumps_structure.job = MUMPSSolver::JOB_INIT;
      this->mumps_structure.comm_fortran = USE_COMM_WORLD;
      dmumps_c(&this->mumps_structure);
      // control parameters
      this->mumps_structure.icntl[0] = -1;
      this->mumps_structure.icntl[1] = -1;
      this->mumps_structure.icntl[2] = -1;
      this->mumps_structure.icntl[3] = 0;
      //this->mumps_structure.icntl[5] = 0; // no scaling
      //this->mumps_structure.icntl[7] = 0; // no scaling

      this->mumps_structure.icntl[12] = 1;
      this->mumps_structure.icntl[23] = 1; // ICNTL(24) controls the detection of “null pivot rows”

      /*
      // debug for MUMPS team
      this->mumps_structure.icntl[1] = 6; // ICNTL(2)=6
      this->mumps_structure.icntl[2] = 6; // ICNTL(3)=6
      this->mumps_structure.icntl[3] = 6; // ICNTL(4)=2
       */

      this->row_indices.reserve(number_nonzeros);
      this->column_indices.reserve(number_nonzeros);
   }

   MUMPSSolver::~MUMPSSolver() {
      this->mumps_structure.job = MUMPSSolver::JOB_END;
      dmumps_c(&this->mumps_structure);
   }
   
   void MUMPSSolver::do_symbolic_analysis(const SymmetricMatrix<size_t, double>& matrix) {
      this->mumps_structure.job = MUMPSSolver::JOB_ANALYSIS;
      this->mumps_structure.n = static_cast<int>(matrix.dimension());
      this->mumps_structure.nnz = static_cast<int>(matrix.number_nonzeros());
      this->mumps_structure.a = nullptr;
      this->save_sparsity_to_local_format(matrix);
      // connect the local sparsity with the pointers in the structure
      this->mumps_structure.irn = this->row_indices.data();
      this->mumps_structure.jcn = this->column_indices.data();
      this->mumps_structure.a = nullptr;
      dmumps_c(&this->mumps_structure);
   }

   void MUMPSSolver::do_numerical_factorization(const SymmetricMatrix<size_t, double>& matrix) {
      this->mumps_structure.job = MUMPSSolver::JOB_FACTORIZATION;
      this->mumps_structure.a = const_cast<double*>(matrix.data_pointer());
      dmumps_c(&this->mumps_structure);
   }

   void MUMPSSolver::solve_indefinite_system(const SymmetricMatrix<size_t, double>& /*matrix*/, const Vector<double>& rhs, Vector<double>& result) {
      result = rhs;
      this->mumps_structure.rhs = result.data();
      this->mumps_structure.job = MUMPSSolver::JOB_SOLVE;
      dmumps_c(&this->mumps_structure);
   }

   std::tuple<size_t, size_t, size_t> MUMPSSolver::get_inertia() const {
      const size_t number_negative_eigenvalues = this->number_negative_eigenvalues();
      const size_t number_zero_eigenvalues = this->number_zero_eigenvalues();
      const size_t number_positive_eigenvalues = static_cast<size_t>(this->mumps_structure.n) - (number_negative_eigenvalues + number_zero_eigenvalues);
      return std::make_tuple(number_positive_eigenvalues, number_negative_eigenvalues, number_zero_eigenvalues);
   }

   size_t MUMPSSolver::number_negative_eigenvalues() const {
      // INFOG(12)
      return static_cast<size_t>(this->mumps_structure.infog[11]);
   }

   size_t MUMPSSolver::number_zero_eigenvalues() const {
      // INFOG(28)
      return static_cast<size_t>(this->mumps_structure.infog[27]);
   }

   bool MUMPSSolver::matrix_is_singular() const {
      return (this->number_zero_eigenvalues() > 0);
   }

   size_t MUMPSSolver::rank() const {
      return this->dimension - this->number_zero_eigenvalues();
   }

   void MUMPSSolver::save_sparsity_to_local_format(const SymmetricMatrix<size_t, double>& matrix) {
      // build the internal matrix representation
      this->row_indices.clear();
      this->column_indices.clear();
      for (const auto [row_index, column_index, _]: matrix) {
         this->row_indices.emplace_back(static_cast<int>(row_index + this->fortran_shift));
         this->column_indices.emplace_back(static_cast<int>(column_index + this->fortran_shift));
      }
   }
} // namespace
