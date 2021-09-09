#include <iostream>
#include <cassert>
#include "PardisoSolver.hpp"

extern "C" {
   void pardisoinit(long int* pt, const int* mtype, const int* solver, int* iparm, double* dparm, int* error);
void pardiso(long int* pt, const int* maxfct, const int* mnum, const int* mtype, const int* phase, const int* n, double* matrix, int* column_start,
      int* row_index, int* perm, const int* nrhs, int* iparm, const int* msglvl, const double* rhs, double* x, int* error, double* dparm);
void pardiso_chkmatrix(const int* mtype, const int* n, double* matrix, int* column_start, int* row_index, int* error);
void pardiso_chkvec(const int* n, const int* nrhs, const double* rhs, int* error);
}

PardisoSolver::PardisoSolver(size_t dimension): LinearSolver(dimension) {
   /* Setup Pardiso control parameters */
   int error = 0;
   pardisoinit(this->pt.data(), &this->mtype, &this->solver, this->iparm.data(), this->dparm.data(), &error);
   assert(error != -10 && "No Pardiso license file found");
   assert(error != -11 && "Pardiso license is expired");
   assert(error != -12 && "Wrong username or hostname");
}

void PardisoSolver::factorize(CSCSymmetricMatrix& matrix) {
   // general factorization method: symbolic factorization and numerical factorization
   this->do_symbolic_factorization(matrix);
   this->do_numerical_factorization(matrix);
}

void PardisoSolver::do_symbolic_factorization(CSCSymmetricMatrix& matrix) {
   // force zero diagonal elements to appear
   matrix.force_explicit_diagonal_elements();

   // convert matrix from 0-based C-notation to Fortran 1-based notation
   for (size_t i = 0; i < this->dimension + 1; i++) {
      matrix.column_start[i]++;
   }
   for (size_t i = 0; i < matrix.number_nonzeros; i++) {
      matrix.row_index[i]++;
   }

   /* check the consistency of the given matrix. Use this functionality only for debugging purposes */
   int error;
   const int dimension = static_cast<int>(this->dimension);
   pardiso_chkmatrix(&this->mtype, &dimension, matrix.matrix.data(), matrix.column_start.data(), matrix.row_index.data(), &error);
   assert(error == 0 && "Consistency error in Pardiso matrix");

   /* Reordering and Symbolic Factorization.  This step also allocates all memory that is necessary for the factorization */
   const int phase = static_cast<int>(ANALYSIS);
   this->iparm[2] = 4; // Number of processors
   pardiso(this->pt.data(), &this->maxfct, &this->mnum, &this->mtype, &phase, &dimension, matrix.matrix.data(),
         matrix.column_start.data(), matrix.row_index.data(), nullptr, &this->nrhs, this->iparm.data(), &this->msglvl, nullptr, nullptr, &error,
         this->dparm.data());
   assert(error == 0 && "Error during Pardiso symbolic factorization");
}

void PardisoSolver::do_numerical_factorization(CSCSymmetricMatrix& matrix) {
   /* Numerical factorization */
   const int phase = static_cast<int>(NUMERICAL_FACTORIZATION);
   this->iparm[32] = 1; /* compute determinant */
   int error;
   const int dimension = static_cast<int>(this->dimension);
   pardiso(this->pt.data(), &this->maxfct, &this->mnum, &this->mtype, &phase, &dimension, matrix.matrix.data(),
         matrix.column_start.data(), matrix.row_index.data(), nullptr, &this->nrhs, this->iparm.data(), &this->msglvl, nullptr, nullptr, &error,
         this->dparm.data());
   assert(error == 0 && "Error occurred during Pardiso numerical factorization");
}

void PardisoSolver::solve(CSCSymmetricMatrix& matrix, const std::vector<double>& rhs, std::vector<double>& result) {
   // check the given vectors for infinite and NaN values
   int error;
   const int dimension = static_cast<int>(this->dimension);
   pardiso_chkvec(&dimension, &this->nrhs, rhs.data(), &error);
   assert(error == 0 && "Error in Pardiso right-hand side");

   /* Back substitution and iterative refinement */
   int phase = static_cast<int>(SOLVE);
   this->iparm[7] = 1; /* Max numbers of iterative refinement steps. */
   pardiso(this->pt.data(), &maxfct, &mnum, &this->mtype, &phase, &dimension, matrix.matrix.data(), matrix.column_start.data(),
         matrix.row_index.data(), nullptr, &nrhs, this->iparm.data(), &msglvl, rhs.data(), result.data(), &error, this->dparm.data());
   assert(error == 0 && "Error during Pardiso solve");

   /* Convert matrix back to 0-based C-notation */
   for (size_t i = 0; i < this->dimension + 1; i++) {
      matrix.column_start[i]--;
   }
   for (size_t i = 0; i < matrix.number_nonzeros; i++) {
      matrix.row_index[i]--;
   }

   /* Release internal memory. */
   phase = static_cast<int>(MEMORY_DEALLOCATION);
   pardiso(this->pt.data(), &this->maxfct, &this->mnum, &this->mtype, &phase, &dimension, nullptr, matrix.column_start.data(),
         matrix.row_index.data(), nullptr, &this->nrhs, this->iparm.data(), &this->msglvl, nullptr, nullptr, &error, this->dparm.data());
}

std::tuple<size_t, size_t, size_t> PardisoSolver::get_inertia() const {
   const size_t rank = this->rank();
   const size_t number_positive_eigenvalues = this->number_positive_eigenvalues();
   const size_t number_negative_eigenvalues = this->number_negative_eigenvalues();
   const size_t number_zero_eigenvalues = this->dimension - rank;
   return std::make_tuple(number_positive_eigenvalues, number_negative_eigenvalues, number_zero_eigenvalues);
}

size_t PardisoSolver::number_positive_eigenvalues() const {
   return static_cast<size_t>(this->iparm[21]);
}

size_t PardisoSolver::number_negative_eigenvalues() const {
   return static_cast<size_t>(this->iparm[22]);
}

bool PardisoSolver::matrix_is_singular() const {
   return (this->rank() == 0);
}

size_t PardisoSolver::rank() const {
   return this->number_positive_eigenvalues() + this->number_negative_eigenvalues();
}