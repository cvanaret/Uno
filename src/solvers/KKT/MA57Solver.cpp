#include <iostream>
#include "MA57Solver.hpp"
#include "Utils.hpp"

extern "C" {
	extern void ma57id_(double cntl[], int icntl[]);
	extern void ma57ad_(int* n, int* ne, int irn[], int jcn[], int* lkeep, int keep[], int iwork[], int icntl[],
		int info[], double rinfo[]);
	extern void ma57bd_(int* n, int* ne, double a[], double fact[], int* lfact, int ifact[], int* lifact, int* lkeep,
		int keep[], int iwork[], int icntl[], double cntl[], int info[], double rinfo[]);
	extern void ma57cd_(int* job, int* n, double fact[], int* lfact, int ifact[], int* lifact, int* nrhs, double rhs[],
		int* lrhs, double work[], int* lwork, int iwork[], int icntl[], int info[]);
}

MA57Solver::MA57Solver(): use_fortran(1), cntl_(5), icntl_(20), info_(40), rinfo_(20) {
}

LocalSolution MA57Solver::solve(COOMatrix& matrix, std::vector<double> rhs) {
	int n = matrix.size;
	int nnz = matrix.number_nonzeros;
	
	/* initialize */
	ma57id_(this->cntl_.data(), this->icntl_.data());
	
	/* analyze sparsity pattern */
	int lkeep = 5*n + nnz + std::max(n, nnz) + 42;
	std::vector<int> keep(lkeep);
	std::vector<int> iwork(5*n);
	std::vector<int> row_indices(matrix.row_indices);
	std::vector<int> column_indices(matrix.column_indices);
	// shift the indices for Fortran
	for (unsigned int k = 0; k < row_indices.size(); k++) {
		row_indices[k] += this->use_fortran;
	}
	for (unsigned int k = 0; k < column_indices.size(); k++) {
		column_indices[k] += this->use_fortran;
	}
	ma57ad_(&n, &nnz, row_indices.data(), column_indices.data(), &lkeep, keep.data(), iwork.data(),
		this->icntl_.data(), this->info_.data(), this->rinfo_.data());

	/* factorize */
	int lfact = 1.2*this->info_[8];
	std::vector<double> fact(lfact);
	int lifact = 1.2*this->info_[9];
	std::vector<int> ifact(lifact);
	ma57bd_(&n, &nnz, matrix.matrix.data(), fact.data(), &lfact, ifact.data(), &lifact, &lkeep, keep.data(),
		iwork.data(), this->icntl_.data(), this->cntl_.data(), this->info_.data(), this->rinfo_.data());;
	
	/* solve */
	int job = 1;
	int nrhs = 1; // number of right hand side being solved
	int lrhs = n; // integer, length of rhs
	int lwork = 1.2*n*nrhs; // length of w; lw>=n*nrhs
	std::vector<double> work(lwork);
	// copy the RHS in the solution
	std::vector<double> x(rhs);
	ma57cd_(&job, &n, fact.data(), &lfact, ifact.data(), &lifact, &nrhs, x.data(), &lrhs, work.data(),
		&lwork, iwork.data(), this->icntl_.data(), this->info_.data());
	
	LocalSolution solution(x, 1, 1);
	return solution;
}

void test() {
	/*
	A[0][0] = 2; A[0][1] = 3;
	A[1][2] = 4; A[1][4] = 6;
	A[2][2] = 1; A[2][3] = 5;
	A[4][4] = 1;
	*/

	int n = 5;
	COOMatrix matrix(n, 0);
	matrix.add_term(2., 0, 0);
	matrix.add_term(3., 0, 1);
	matrix.add_term(4., 1, 2);
	matrix.add_term(6., 1, 4);
	matrix.add_term(1., 2, 2);
	matrix.add_term(5., 2, 3);
	matrix.add_term(1., 4, 4);

	/* right hand side */
	std::vector<double> rhs = {8., 45., 31., 15., 17.};

	MA57Solver s;
	LocalSolution solution = s.solve(matrix, rhs);

	std::cout << "solution";
	for (unsigned int k = 0; k < solution.x.size(); k++) {
		std::cout << " " << solution.x[k];
	}
	std::cout << "\n";
	return;
}
