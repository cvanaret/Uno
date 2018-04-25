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

MA57Solver::MA57Solver() {
	this->lkeep_ = 100;
	this->lrhs_ = 10;
	this->lwork_ = 12;
	this->lfact_ = 1000;
	this->lifact_ = 500;
}

void hessian_sparsity_from_ampl(std::vector<int>& hessian_column_start, std::vector<int>& hessian_row_number) {
	int n = hessian_column_start.size() - 1;
	this->hessian_nnz = hessian_row_number.size();
	
	/* MA57 sparsity pattern */
	this->hessian_row_number = std::vector<int>(this->hessian_nnz);
	this->hessian_column_number = std::vector<int>(this->hessian_nnz);
	
	int use_fortran = 1;
	for (int j = 0; j < n; j++) {
		for (int k = matrix_column_start[j]-use_fortran; k < matrix_column_start[j+1]-use_fortran; k++) {
			int i = matrix_row_number[k];
			
			this->hessian_row_number.push_back(i);
			this->hessian_column_number.push_back(j + use_fortran);
		}
	}
	return;
}

void test_ma57() {
	std::vector<int> ampl_hessian_column_start = {1, 2, 4, 7, 11, 11};
	std::vector<int> ampl_hessian_row_number = {1, 1, 2, 1, 2, 3, 1, 2, 3, 4};
	hessian_sparsity_from_ampl(ampl_hessian_column_start, ampl_hessian_row_number);
	/*
	Rows:    1 1 2 1 2 3 1 2 3 4
	Columns: 1 2 2 3 3 3 4 4 4 4
	*/
}
