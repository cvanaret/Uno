#include "MUMPSSolver.hpp"
#include "mpi.h"

#define USE_COMM_WORLD -987654

MUMPSSolver::MUMPSSolver(size_t dimension, size_t number_nonzeros): SymmetricIndefiniteLinearSolver<size_t, double>(dimension),
      COO_matrix(dimension, number_nonzeros, false) {
   this->mumps_structure.sym = MUMPSSolver::GENERAL_SYMMETRIC;
   this->mumps_structure.par = 1;
   this->mumps_structure.job = MUMPSSolver::JOB_INIT;
   this->mumps_structure.comm_fortran = USE_COMM_WORLD;
   dmumps_c(&this->mumps_structure);
   // control parameters
   this->mumps_structure.icntl[0] = -1;
   this->mumps_structure.icntl[1] = -1;
   this->mumps_structure.icntl[2] = -1;
   this->mumps_structure.icntl[3] = 0;
}

MUMPSSolver::~MUMPSSolver() {
   this->mumps_structure.job = MUMPSSolver::JOB_END;
   dmumps_c(&this->mumps_structure);
}

void MUMPSSolver::factorize(const SymmetricMatrix<size_t, double>& matrix) {
   // general factorization method: symbolic factorization and numerical factorization
   this->do_symbolic_factorization(matrix);
   this->do_numerical_factorization(matrix);
}

void MUMPSSolver::do_symbolic_factorization(const SymmetricMatrix<size_t, double>& matrix) {
   this->save_matrix_to_local_format(matrix);
   this->mumps_structure.n = static_cast<int>(matrix.dimension);
   this->mumps_structure.nnz = static_cast<int>(matrix.number_nonzeros);
   this->mumps_structure.job = MUMPSSolver::JOB_ANALYSIS;
   // connect the local COO matrix with the pointers in the structure
   this->mumps_structure.irn = this->COO_matrix.row_indices_pointer();
   this->mumps_structure.jcn = this->COO_matrix.column_indices_pointer();
   dmumps_c(&this->mumps_structure);
}

void MUMPSSolver::do_numerical_factorization(const SymmetricMatrix<size_t, double>& /*matrix*/) {
   this->mumps_structure.job = MUMPSSolver::JOB_FACTORIZATION;
   this->mumps_structure.a = this->COO_matrix.data_pointer();
   dmumps_c(&this->mumps_structure);
}

void MUMPSSolver::solve_indefinite_system(const SymmetricMatrix<size_t, double>& /*matrix*/, const Vector<double>& rhs, Vector<double>& result) {
   result = rhs;
   this->mumps_structure.rhs = result.data();
   this->mumps_structure.job = MUMPSSolver::JOB_SOLVE;
   dmumps_c(&this->mumps_structure);
}

std::tuple <size_t, size_t, size_t> MUMPSSolver::get_inertia() const {
   // TODO
   return {0, 0, 0};
}

size_t MUMPSSolver::number_negative_eigenvalues() const {
   return static_cast<size_t>(this->mumps_structure.infog[11]);
}

bool MUMPSSolver::matrix_is_singular() const {
   // TODO
   return false;
}

size_t MUMPSSolver::rank() const {
   // TODO
   return 0;
}

void MUMPSSolver::save_matrix_to_local_format(const SymmetricMatrix<size_t, double>& matrix) {
   // build the internal matrix representation
   this->COO_matrix.reset();
   for (const auto [row_index, column_index, element]: matrix) {
      this->COO_matrix.insert(element, static_cast<int>(row_index + this->fortran_shift), static_cast<int>(column_index + this->fortran_shift));
   }
}
