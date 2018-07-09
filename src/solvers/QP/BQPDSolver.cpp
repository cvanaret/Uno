#include <cmath>
#include "BQPDSolver.hpp"

#define BIG 1e30

extern "C" {
	/* fortran common block used in bqpd/bqpd.f */
	extern struct {
		int kk, ll, kkk, lll, mxws, mxlws;
	} wsc_;

	/* fortran common for inertia correction in wdotd */
	extern struct {
		double alpha;
	} kktalphac_;

	extern void bqpd_(int *n, int *m, int *k, int *kmax, double *a, int *la, double *x, double *bl, double *bu, 
			   double *f, double *fmin, double *g, double *r, double *w, double *e, int *ls, double *alp, int *lp,
			   int *mlp, int *peq, double *ws, int *lws, int *mode, int *ifail, int *info, int *iprint, int *nout);
}

BQPDSolver::BQPDSolver(std::vector<int>& hessian_column_start, std::vector<int>& hessian_row_number): QPSolver(), LPSolver(),
		kmax_(500), mlp_(1000), mxwk0_(2000000), mxiwk0_(500000), info_(100), alp_(mlp_), lp_(mlp_),
		hessian_column_start(hessian_column_start), hessian_row_number(hessian_row_number) {
	/* preallocate a bunch of stuff */
	this->k_ = 0;
	/* warm start: mode = 2 */
	this->mode_ = 0;
	
	this->iprint_ = 0;
	this->nout_ = 6;
	
	this->fmin_ = -1e20;
	
	this->hessian_nnz_ = hessian_row_number.size();
	this->nhr_ = this->hessian_nnz_;
	
	kktalphac_.alpha = 0;
	this->use_fortran = 1;
}

void BQPDSolver::allocate(int n, int m) {
	this->n_ = n;
	this->m_ = m;
	this->ls_.resize(n + m);
	this->w_.resize(n + m);
	this->gradient_solution_.resize(n);
	this->residuals_.resize(n + m);
	this->e_.resize(n + m);
	this->nhi_ = hessian_nnz_ + n + 3;
	this->mxws_ = nhr_ + kmax_*(kmax_+9)/2 + 2*n + m + mxwk0_;
	this->mxlws_ = nhi_ + kmax_ + mxiwk0_;
	this->ws_.resize(mxws_);
	this->lws_.resize(mxlws_);
	
	/* active set */
	for (int i = 0; i < n+m; i++) {
		this->ls_[i] = i+1;
	}
	
	/* Hessian sparsity */
	int last_value = this->hessian_column_start[this->hessian_column_start.size()-1];
	for (int j = this->hessian_column_start.size()-1; j < n; j++) {
		this->hessian_column_start.push_back(last_value);
	}
	
	this->lws_[0] = this->hessian_nnz_ + 1;
	for(int i = 0; i < this->hessian_nnz_; i++) {
		this->lws_[i+1] = hessian_row_number[i] + this->use_fortran;
	}
	for(int i = 0; i < this->n_ + 1; i++) {
		this->lws_[this->hessian_nnz_ + i + 1] = hessian_column_start[i] + this->use_fortran;
	}
	
	/* initialize wsc_ common block (Hessian & workspace for bqpd) */
	wsc_.kk = this->nhr_;
	wsc_.ll = this->nhi_;
	wsc_.mxws = this->mxws_;
	wsc_.mxlws = this->mxlws_;
	return;
}

Status int_to_status(int ifail) {
	if (ifail < 0 || 10 <= ifail) {
		throw std::length_error("BQPDSolver.int_to_status: ifail does not belong to [0, 9]");
	}
	Status status = static_cast<Status>(ifail);
	return status;
}

LocalSolution BQPDSolver::solve(QP& qp, std::vector<double>& x) {
	/* Hessian */
	for (unsigned int i = 0; i < qp.hessian.number_nonzeros; i++) {
		this->ws_[i] = qp.hessian.matrix[i];
	}
	
	/* Jacobian */
	// TODO preallocate
	std::vector<double> jacobian;
	std::vector<int> jacobian_sparsity;
	jacobian_sparsity.push_back(0); // header-related, will be modified later on
	build_jacobian(jacobian, jacobian_sparsity, qp.objective);
	for (int j = 0; j < qp.number_constraints; j++) {
		build_jacobian(jacobian, jacobian_sparsity, qp.constraints[j]);
	}
	/* Jacobian header */
	jacobian_sparsity[0] = jacobian_sparsity.size();
	unsigned int total_size = 1;
	jacobian_sparsity.push_back(total_size);
	total_size += qp.objective.size();
	jacobian_sparsity.push_back(total_size);
	for (int j = 0; j < qp.number_constraints; j++) {
		total_size += qp.constraints[j].size();
		jacobian_sparsity.push_back(total_size);
	}
	
	/* bounds */
	std::vector<double> lb(qp.number_variables + qp.number_constraints);
	std::vector<double> ub(qp.number_variables + qp.number_constraints);
	for (int i = 0; i < qp.number_variables; i++) {
		lb[i] = qp.variable_lb[i];
		ub[i] = qp.variable_ub[i];
		if (lb[i] == -INFINITY) {
			lb[i] = -BIG;
		}
		if (ub[i] == INFINITY) {
			ub[i] = BIG;
		}
	}
	for (int j = 0; j < qp.number_constraints; j++) {
		lb[qp.number_variables + j] = qp.constraint_lb[j];
		ub[qp.number_variables + j] = qp.constraint_ub[j];
	}
	
	/* call BQPD */
	bqpd_(&this->n_, &this->m_, &this->k_, &this->kmax_, jacobian.data(), jacobian_sparsity.data(), x.data(),
		lb.data(), ub.data(), &this->f_solution_, &this->fmin_, this->gradient_solution_.data(),
		this->residuals_.data(), this->w_.data(), this->e_.data(), this->ls_.data(), this->alp_.data(),
		this->lp_.data(), &this->mlp_, &this->peq_solution_, this->ws_.data(), this->lws_.data(), &this->mode_,
		&this->ifail_, this->info_.data(), &this->iprint_, &this->nout_);

	/* project solution into bounds: it's a ray! */
	for (unsigned int i = 0; i < x.size(); i++) {
		if (x[i] < qp.variable_lb[i]) {
			x[i] = qp.variable_lb[i];
		}
		else if (qp.variable_ub[i] < x[i]) {
			x[i] = qp.variable_ub[i];
		}
	}

	LocalSolution solution = this->generate_solution(x);
	return solution;
}

LocalSolution BQPDSolver::solve(LP& lp, std::vector<double>& x) {
	int kmax = 0;
	
	/* Jacobian */
	// TODO preallocate
	std::vector<double> jacobian;
	std::vector<int> jacobian_sparsity;
	jacobian_sparsity.push_back(0); // header-related, to be modified later on
	build_jacobian(jacobian, jacobian_sparsity, lp.objective);
	for (int j = 0; j < lp.number_constraints; j++) {
		build_jacobian(jacobian, jacobian_sparsity, lp.constraints[j]);
	}
	/* Jacobian header */
	jacobian_sparsity[0] = jacobian_sparsity.size();
	unsigned int total_size = 1;
	jacobian_sparsity.push_back(total_size);
	total_size += lp.objective.size();
	jacobian_sparsity.push_back(total_size);
	for (int j = 0; j < lp.number_constraints; j++) {
		total_size += lp.constraints[j].size();
		jacobian_sparsity.push_back(total_size);
	}

	/* bounds */
	std::vector<double> lb(lp.number_variables + lp.number_constraints);
	std::vector<double> ub(lp.number_variables + lp.number_constraints);
	for (int i = 0; i < lp.number_variables; i++) {
		lb[i] = lp.variable_lb[i];
		ub[i] = lp.variable_ub[i];
		if (lb[i] == -INFINITY) {
			lb[i] = -BIG;
		}
		if (ub[i] == INFINITY) {
			ub[i] = BIG;
		}
	}
	for (int j = 0; j < lp.number_constraints; j++) {
		lb[lp.number_variables + j] = lp.constraint_lb[j];
		ub[lp.number_variables + j] = lp.constraint_ub[j];
	}

	/* call BQPD */
	bqpd_(&this->n_, &this->m_, &this->k_, &kmax, jacobian.data(), jacobian_sparsity.data(), x.data(),
		lb.data(), ub.data(), &this->f_solution_, &this->fmin_, this->gradient_solution_.data(),
		this->residuals_.data(), this->w_.data(), this->e_.data(), this->ls_.data(), this->alp_.data(),
		this->lp_.data(), &this->mlp_, &this->peq_solution_, this->ws_.data(), this->lws_.data(), &this->mode_,
		&this->ifail_, this->info_.data(), &this->iprint_, &this->nout_);

	/* project solution into bounds: it's a ray! */
	for (unsigned int i = 0; i < x.size(); i++) {
		if (x[i] < lp.variable_lb[i]) {
			x[i] = lp.variable_lb[i];
		}
		else if (lp.variable_ub[i] < x[i]) {
			x[i] = lp.variable_ub[i];
		}
	}

	LocalSolution solution = this->generate_solution(x);
	return solution;
}

LocalSolution BQPDSolver::generate_solution(std::vector<double>& x) {
	/* generate the result */
	LocalSolution solution(x, this->n_, this->m_);
	solution.status = int_to_status(this->ifail_);
	solution.norm = norm_inf(x);
	solution.objective = this->f_solution_;
	std::vector<ConstraintFeasibility> status(this->m_);
	
	/* active constraints */
	for (int j = 0; j < this->n_-this->k_; j++) {
		int index = std::abs(this->ls_[j])-1;
		
		if (this->ls_[j] < 0) { /* upper bound active */
			solution.multipliers[index] = -this->residuals_[index];
			solution.active_set.at_upper_bound.push_back(index);
		}
		else { /* lower bound active */
			solution.multipliers[index] = this->residuals_[index];
			solution.active_set.at_lower_bound.push_back(index);
		}
		if (this->n_ <= index) {
			solution.constraint_partition.feasible_set.push_back(index - this->n_);
			status[index - this->n_] = FEASIBLE;
		}
	}
	
	/* inactive constraints */
	for (int j = this->n_-this->k_; j < this->n_+this->m_; j++) {
		int index = std::abs(this->ls_[j])-1;
		solution.multipliers[index] = 0.;
		
		if (this->n_ <= index) { // general constraints
			if (this->residuals_[index] < 0.) { // infeasible constraint
				solution.constraint_partition.infeasible_set.push_back(index - this->n_);
				if (this->ls_[j] < 0) { // upper bound violated
					status[index - this->n_] = INFEASIBLE_UPPER;
				}
				else { // lower bound violated
					status[index - this->n_] = INFEASIBLE_LOWER;
				}
			}
			else { // feasible constraint
				solution.constraint_partition.feasible_set.push_back(index - this->n_);
				status[index - this->n_] = FEASIBLE;
			}
		}
	}
	solution.constraint_partition.status = status;

	return solution;
}

void BQPDSolver::build_jacobian(std::vector<double>& full_jacobian, std::vector<int>& full_jacobian_sparsity, std::map<int,double>& jacobian) {
	for (std::map<int,double>::iterator it = jacobian.begin(); it != jacobian.end(); it++) {
		int variable_index = it->first;
		double derivative = it->second;
		
		full_jacobian.push_back(derivative);
		full_jacobian_sparsity.push_back(variable_index + this->use_fortran);
	}
	return;
}
