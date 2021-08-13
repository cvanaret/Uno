#include <iostream>
#include <cassert>
#include "MA57Solver.hpp"
#include "linear_algebra/Vector.hpp"

extern "C" {
// MA57
// default values of controlling parameters
void ma57id_(double cntl[], int icntl[]);
// symbolic factorization
void
ma57ad_(const int* n, const int* ne, const int irn[], const int jcn[], const int* lkeep, int keep[], int iwork[], int icntl[], int info[], double
rinfo[]);
// numerical factorization
void ma57bd_(int* n, int* ne, const double a[], /* out */ double fact[], const int* lfact, /* out */ int ifact[], const int* lifact,
      const int* lkeep, const int keep[], int iwork[], int icntl[], double cntl[], /* out */ int info[], /* out */ double rinfo[]);
// linear system solve without iterative refinement
void ma57cd_(const int* job, int* n, double fact[], int* lfact, int ifact[], int* lifact, const int* nrhs, double rhs[], int* lrhs, double work[],
      int* lwork, int iwork[], int icntl[], int info[]);
// linear system solve with iterative refinement
void ma57dd_(const int* job, int* n, int* ne, const double a[], const int irn[], const int jcn[], double fact[], int* lfact, int ifact[], int* lifact,
      const double rhs[], double x[], double resid[], double work[], int iwork[], int icntl[],
      double cntl[], int info[], double rinfo[]);
}

MA57Solver::MA57Solver() : LinearSolver() {
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
   const int n = matrix.dimension;
   const int nnz = (int) matrix.number_nonzeros;

   /* sparsity pattern */
   const int lkeep = 5 * n + nnz + std::max(n, nnz) + 42;
   std::vector<int> keep(lkeep);
   std::vector<int> iwork(5 * n);

   /* reindex the matrix (Fortran compliance) */
   for (int k = 0; k < matrix.number_nonzeros; k++) {
      matrix.row_indices[k]++;
      matrix.column_indices[k]++;
   }

   /* symbolic factorization */
   std::array<int, 40> info{};
   ma57ad_(/* const */ &n,
         /* const */ &nnz,
         /* const */ matrix.row_indices.data(),
         /* const */ matrix.column_indices.data(),
         /* const */ &lkeep,
         /* const */ keep.data(),
         /* out */ iwork.data(),
         /* const */ this->icntl.data(),
         /* out */ info.data(),
         /* out */ this->rinfo.data());

   /* reindex the matrix (Fortran compliance) */
   for (int k = 0; k < matrix.number_nonzeros; k++) {
      matrix.row_indices[k]--;
      matrix.column_indices[k]--;
   }

   // TODO check that info[0] == 0
   int lfact = 2 * info[8];
   std::vector<double> fact(lfact);
   int lifact = 2 * info[9];
   std::vector<int> ifact(lifact);
   // build the factorization object
   this->factorization = {n, nnz, fact, lfact, ifact, lifact, lkeep, keep, iwork, info};
}

void MA57Solver::do_numerical_factorization(COOSymmetricMatrix& matrix) {
   assert(this->factorization.n == matrix.dimension && "MA57Solver: the dimensions do not match");
   assert(this->factorization.nnz == matrix.number_nonzeros && "MA57Solver: the numbers of nonzeros do not match");

   /* numerical factorization */
   ma57bd_(&this->factorization.n,
         (int*) &this->factorization.nnz,
         /* const */ matrix.matrix.data(),
         /* out */ this->factorization.fact.data(),
         /* const */ &this->factorization.lfact,
         /* out*/ this->factorization.ifact.data(),
         /* const */ &this->factorization.lifact,
         /* const */ &this->factorization.lkeep,
         /* const */ this->factorization.keep.data(), this->factorization.iwork.data(), this->icntl.data(), this->cntl.data(),
         /* out */ this->factorization.info.data(),
         /* out */ this->rinfo.data());
}

std::vector<double> MA57Solver::solve(COOSymmetricMatrix& matrix, const std::vector<double>& rhs) {
   /* solve */
   int n = this->factorization.n;
   int lrhs = n; // integer, length of rhs
   int lwork = 1.2 * n * this->nrhs; // length of w; lw >= n*nrhs
   std::vector<double> work(lwork);

   // allocate a vector for the solution
   std::vector<double> x(n);

   // solve the linear system
   if (this->use_iterative_refinement) {
      std::vector<double> residuals(n);
      assert(false && "TODO in MA57Solver::solve: reindex matrix");

      ma57dd_(&this->job, &n, (int*) &this->factorization.nnz, matrix.matrix.data(), matrix.row_indices.data(), matrix.column_indices.data(),
            this->factorization.fact.data(), &this->factorization.lfact, this->factorization.ifact.data(), &this->factorization.lifact,
            rhs.data(), x.data(), residuals.data(), work.data(), this->factorization.iwork.data(), this->icntl.data(),
            this->cntl.data(), this->factorization.info.data(), this->rinfo.data());
   }
   else {
      // copy rhs into x
      copy_from(x, rhs);

      ma57cd_(&this->job, &n, this->factorization.fact.data(), &this->factorization.lfact, this->factorization.ifact.data(),
            &this->factorization.lifact, &this->nrhs, x.data(), &lrhs, work.data(), &lwork, this->factorization.iwork.data(),
            this->icntl.data(), this->factorization.info.data());
   }
   return x;
}

std::tuple<int, int, int> MA57Solver::get_inertia() const {
   const int rank = this->rank();
   const int number_negative_eigenvalues = (int) this->number_negative_eigenvalues();
   const int number_positive_eigenvalues = rank - number_negative_eigenvalues;
   const int number_zero_eigenvalues = (int) this->factorization.n - rank;
   return std::make_tuple(number_positive_eigenvalues, number_negative_eigenvalues, number_zero_eigenvalues);
}

size_t MA57Solver::number_negative_eigenvalues() const {
   return this->factorization.info[23];
}

bool MA57Solver::matrix_is_singular() const {
   return (this->factorization.info[0] == 4);
}

int MA57Solver::rank() const {
   return this->factorization.info[24];
}