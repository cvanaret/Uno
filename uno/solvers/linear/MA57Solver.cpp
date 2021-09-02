#include <iostream>
#include <cassert>
#include "MA57Solver.hpp"
#include "linear_algebra/Vector.hpp"

extern "C" {
// MA57
// default values of controlling parameters
void ma57id_(double cntl[], int icntl[]);
// symbolic factorization
void ma57ad_(const int* n, const int* ne, const int irn[], const int jcn[], const int* lkeep, int keep[], int iwork[], int icntl[], int info[], double
rinfo[]);
// numerical factorization
void ma57bd_(const int* n, int* ne, const double a[], /* out */ double fact[], const int* lfact, /* out */ int ifact[], const int* lifact,
      const int* lkeep, const int keep[], int iwork[], int icntl[], double cntl[], /* out */ int info[], /* out */ double rinfo[]);
// linear system solve without iterative refinement
void ma57cd_(const int* job, const int* n, double fact[], int* lfact, int ifact[], int* lifact, const int* nrhs, double rhs[], const int* lrhs, double
work[], int* lwork, int iwork[], int icntl[], int info[]);
// linear system solve with iterative refinement
void ma57dd_(const int* job, const int* n, int* ne, const double a[], const int irn[], const int jcn[], double fact[], int* lfact, int ifact[], int*
lifact, const double rhs[], double x[], double resid[], double work[], int iwork[], int icntl[],
      double cntl[], int info[], double rinfo[]);
}

MA57Solver::MA57Solver(size_t dimension) : LinearSolver(), dimension(dimension), iwork(5*dimension),
   lwork((int) (1.2 * dimension)), work(this->lwork), residuals(dimension) {
   /* set the default values of the controlling parameters */
   ma57id_(this->cntl.data(), this->icntl.data());
   // suppress warning messages
   this->icntl[4] = 0;
   // iterative refinement enabled
   this->icntl[8] = 1;
}

void MA57Solver::factorize(COOSymmetricMatrix& matrix) {
   // general factorization method: symbolic factorization and numerical factorization
   this->do_symbolic_factorization(matrix);
   this->do_numerical_factorization(matrix);
}

void MA57Solver::do_symbolic_factorization(COOSymmetricMatrix& matrix) {
   assert(matrix.dimension == this->dimension && "MA57Solver: the dimension of the matrix is inconsistent");
   const int n = (int) matrix.dimension;
   const int nnz = (int) matrix.number_nonzeros;

   /* sparsity pattern */
   const int lkeep = 5 * n + nnz + std::max(n, nnz) + 42;
   std::vector<int> keep(lkeep);

   /* reindex the matrix (Fortran compliance) */
   for (size_t k = 0; k < matrix.number_nonzeros; k++) {
      matrix.row_indices[k]++;
      matrix.column_indices[k]++;
   }

   /* symbolic factorization */
   ma57ad_(/* const */ &n,
         /* const */ &nnz,
         /* const */ matrix.row_indices.data(),
         /* const */ matrix.column_indices.data(),
         /* const */ &lkeep,
         /* const */ keep.data(),
         /* out */ this->iwork.data(),
         /* const */ this->icntl.data(),
         /* out */ this->info.data(),
         /* out */ this->rinfo.data());

   assert(info[0] == 0 && "MA57: the symbolic factorization failed");
   int lfact = 2 * this->info[8];
   std::vector<double> fact(lfact);
   int lifact = 2 * this->info[9];
   std::vector<int> ifact(lifact);

   /* reindex the matrix (Fortran compliance) */
   for (size_t k = 0; k < matrix.number_nonzeros; k++) {
      matrix.row_indices[k]--;
      matrix.column_indices[k]--;
   }

   // build the factorization object
   this->factorization = {nnz, std::move(fact), lfact, std::move(ifact), lifact, lkeep, std::move(keep)};
}

void MA57Solver::do_numerical_factorization(COOSymmetricMatrix& matrix) {
   assert(this->dimension == matrix.dimension && "MA57Solver: the dimension of the matrix is inconsistent");
   assert(this->factorization.nnz == (int) matrix.number_nonzeros && "MA57Solver: the numbers of nonzeros do not match");

   const int n = this->dimension;
   /* numerical factorization */
   ma57bd_(&n,
         (int*) &this->factorization.nnz,
         /* const */ matrix.matrix.data(),
         /* out */ this->factorization.fact.data(),
         /* const */ &this->factorization.lfact,
         /* out*/ this->factorization.ifact.data(),
         /* const */ &this->factorization.lifact,
         /* const */ &this->factorization.lkeep,
         /* const */ this->factorization.keep.data(), this->iwork.data(), this->icntl.data(), this->cntl.data(),
         /* out */ this->info.data(),
         /* out */ this->rinfo.data());
}

void MA57Solver::solve(COOSymmetricMatrix& matrix, const std::vector<double>& rhs, std::vector<double>& result) {
   /* solve */
   const int n = this->dimension;
   const int lrhs = n; // integer, length of rhs

   // solve the linear system
   if (this->use_iterative_refinement) {
      assert(false && "TODO in MA57Solver::solve: reindex matrix");

      ma57dd_(&this->job, &n, (int*) &this->factorization.nnz, matrix.matrix.data(), matrix.row_indices.data(), matrix.column_indices.data(),
            this->factorization.fact.data(), &this->factorization.lfact, this->factorization.ifact.data(), &this->factorization.lifact,
            rhs.data(), result.data(), this->residuals.data(), this->work.data(), this->iwork.data(), this->icntl.data(),
            this->cntl.data(), this->info.data(), this->rinfo.data());
   }
   else {
      // copy rhs into result (overwritten by MA57)
      copy_from(result, rhs);

      ma57cd_(&this->job, &n, this->factorization.fact.data(), &this->factorization.lfact, this->factorization.ifact.data(),
            &this->factorization.lifact, &this->nrhs, result.data(), &lrhs, this->work.data(), &this->lwork, this->iwork.data(),
            this->icntl.data(), this->info.data());
   }
}

std::tuple<int, int, int> MA57Solver::get_inertia() const {
   const int rank = this->rank();
   const int number_negative_eigenvalues = (int) this->number_negative_eigenvalues();
   const int number_positive_eigenvalues = rank - number_negative_eigenvalues;
   const int number_zero_eigenvalues = (int) this->dimension - rank;
   return std::make_tuple(number_positive_eigenvalues, number_negative_eigenvalues, number_zero_eigenvalues);
}

size_t MA57Solver::number_negative_eigenvalues() const {
   return this->info[23];
}

bool MA57Solver::matrix_is_singular() const {
   return (this->info[0] == 4);
}

int MA57Solver::rank() const {
   return this->info[24];
}