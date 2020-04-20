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

/* preallocate a bunch of stuff */
BQPDSolver::BQPDSolver(std::vector<int>& hessian_column_start, std::vector<int>& hessian_row_number) : QPSolver(),
use_fortran(1), kmax_(500), mlp_(1000), mxwk0_(2000000), mxiwk0_(500000), info_(100), alp_(mlp_), lp_(mlp_),
hessian_nnz_(hessian_row_number.size()), nhr_(hessian_row_number.size()), k_(0), mode_(COLD_START),
iprint_(0), nout_(6), fmin_(-1e20),
hessian_column_start(hessian_column_start), hessian_row_number(hessian_row_number) {
    kktalphac_.alpha = 0;
}

void BQPDSolver::allocate(int n, int m) {
    this->n_ = n;
    this->m_ = m;
    this->ls_.resize(n + m); // active set
    this->w_.resize(n + m);
    this->gradient_solution_.resize(n);
    this->residuals_.resize(n + m); // multipliers
    this->e_.resize(n + m);
    this->nhi_ = hessian_nnz_ + n + 3;
    this->mxws_ = nhr_ + kmax_ * (kmax_ + 9) / 2 + 2 * n + m + mxwk0_;
    this->mxlws_ = nhi_ + kmax_ + mxiwk0_;
    this->ws_.resize(mxws_);
    this->lws_.resize(mxlws_);

    /* active set */
    for (int i = 0; i < n + m; i++) {
        this->ls_[i] = i + 1;
    }

    /* Hessian sparsity */
    int last_value = this->hessian_column_start[this->hessian_column_start.size() - 1];
    for (int j = this->hessian_column_start.size() - 1; j < n; j++) {
        this->hessian_column_start.push_back(last_value);
    }

    this->lws_[0] = this->hessian_nnz_ + 1;
    for (int i = 0; i < this->hessian_nnz_; i++) {
        this->lws_[i + 1] = hessian_row_number[i] + this->use_fortran;
    }
    for (int i = 0; i < this->n_ + 1; i++) {
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
    Status status = static_cast<Status> (ifail);
    return status;
}

SubproblemSolution BQPDSolver::solve_QP(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, std::map<int, double>& linear_objective, std::vector<std::map<int, double> >& constraints_jacobian, CSCMatrix& hessian, std::vector<double>& x) {
    /* Hessian */
    for (int i = 0; i < hessian.number_nonzeros(); i++) {
        this->ws_[i] = hessian.matrix[i];
    }
    return this->solve_subproblem(variables_bounds, constraints_bounds, linear_objective, constraints_jacobian, x, this->kmax_);
}

SubproblemSolution BQPDSolver::solve_LP(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, std::map<int, double>& linear_objective, std::vector<std::map<int, double> >& constraints_jacobian, std::vector<double>& x) {
    return this->solve_subproblem(variables_bounds, constraints_bounds, linear_objective, constraints_jacobian, x, 0);
}

SubproblemSolution BQPDSolver::solve_subproblem(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, std::map<int, double>& linear_objective, std::vector<std::map<int, double> >& constraints_jacobian, std::vector<double>& x, int kmax) {
    /* Jacobian */
    // TODO preallocate
    std::vector<double> jacobian;
    std::vector<int> jacobian_sparsity;
    jacobian_sparsity.push_back(0); // header-related, to be modified later on
    build_jacobian(jacobian, jacobian_sparsity, linear_objective);
    for (unsigned int j = 0; j < constraints_jacobian.size(); j++) {
        build_jacobian(jacobian, jacobian_sparsity, constraints_jacobian[j]);
    }
    /* Jacobian header */
    jacobian_sparsity[0] = jacobian_sparsity.size();
    unsigned int total_size = 1;
    jacobian_sparsity.push_back(total_size);
    total_size += linear_objective.size();
    jacobian_sparsity.push_back(total_size);
    for (unsigned int j = 0; j < constraints_jacobian.size(); j++) {
        total_size += constraints_jacobian[j].size();
        jacobian_sparsity.push_back(total_size);
    }

    /* bounds */
    std::vector<double> lb(variables_bounds.size() + constraints_bounds.size());
    std::vector<double> ub(variables_bounds.size() + constraints_bounds.size());
    for (unsigned int i = 0; i < variables_bounds.size(); i++) {
        lb[i] = (variables_bounds[i].lb == -INFINITY) ? -BIG : variables_bounds[i].lb;
        ub[i] = (variables_bounds[i].ub == INFINITY) ? BIG : variables_bounds[i].ub;
    }
    for (unsigned int j = 0; j < constraints_bounds.size(); j++) {
        lb[variables_bounds.size() + j] = constraints_bounds[j].lb;
        ub[variables_bounds.size() + j] = constraints_bounds[j].ub;
    }

    /* call BQPD */
    int mode = (int) this->mode_;
    bqpd_(&this->n_, &this->m_, &this->k_, &kmax, jacobian.data(), jacobian_sparsity.data(), x.data(),
            lb.data(), ub.data(), &this->f_solution_, &this->fmin_, this->gradient_solution_.data(),
            this->residuals_.data(), this->w_.data(), this->e_.data(), this->ls_.data(), this->alp_.data(),
            this->lp_.data(), &this->mlp_, &this->peq_solution_, this->ws_.data(), this->lws_.data(), &mode,
            &this->ifail_, this->info_.data(), &this->iprint_, &this->nout_);

    /* project solution into bounds: it's a ray! */
    for (unsigned int i = 0; i < x.size(); i++) {
        if (x[i] < variables_bounds[i].lb) {
            x[i] = variables_bounds[i].lb;
        }
        else if (variables_bounds[i].ub < x[i]) {
            x[i] = variables_bounds[i].ub;
        }
    }

    SubproblemSolution solution = this->generate_solution(x);
    return solution;
}

SubproblemSolution BQPDSolver::generate_solution(std::vector<double>& x) {
    Multipliers multipliers(this->n_, this->m_);
    SubproblemSolution solution(x, multipliers);

    /* active constraints */
    for (int j = 0; j < this->n_ - this->k_; j++) {
        int index = std::abs(this->ls_[j]) - this->use_fortran;

        if (this->ls_[j] < 0) { /* upper bound active */
            solution.active_set.at_upper_bound.insert(index);
        }
        else { /* lower bound active */
            solution.active_set.at_lower_bound.insert(index);
        }

        if (index < this->n_) {
            if (this->ls_[j] < 0) { /* upper bound active */
                solution.multipliers.upper_bounds[index] = -this->residuals_[index];
            }
            else { /* lower bound active */
                solution.multipliers.lower_bounds[index] = this->residuals_[index];
            }
        }
        else {
            int constraint_index = index - this->n_;
            solution.constraint_partition.feasible.insert(constraint_index);
            solution.constraint_partition.constraint_feasibility[constraint_index] = FEASIBLE;
            solution.multipliers.constraints[constraint_index] = (this->ls_[j] < 0) ? - this->residuals_[index] : this->residuals_[index];
        }
    }

    /* inactive constraints */
    for (int j = this->n_ - this->k_; j < this->n_ + this->m_; j++) {
        int index = std::abs(this->ls_[j]) - this->use_fortran;

        if (this->n_ <= index) { // general constraints
            int constraint_index = index - this->n_;
            if (this->residuals_[index] < 0.) { // infeasible constraint
                solution.constraint_partition.infeasible.insert(constraint_index);
                if (this->ls_[j] < 0) { // upper bound violated
                    solution.constraint_partition.constraint_feasibility[constraint_index] = INFEASIBLE_UPPER;
                }
                else { // lower bound violated
                    solution.constraint_partition.constraint_feasibility[constraint_index] = INFEASIBLE_LOWER;
                }
            }
            else { // feasible constraint
                solution.constraint_partition.feasible.insert(constraint_index);
                solution.constraint_partition.constraint_feasibility[constraint_index] = FEASIBLE;
            }
        }
    }
    solution.status = int_to_status(this->ifail_);
    // phase
    // phase_1_required
    solution.norm = norm_inf(x);
    solution.objective = this->f_solution_;
    return solution;
}

void BQPDSolver::build_jacobian(std::vector<double>& full_jacobian, std::vector<int>& full_jacobian_sparsity, std::map<int, double>& jacobian) {
    for (std::map<int, double>::iterator it = jacobian.begin(); it != jacobian.end(); it++) {
        int i = it->first;
        double derivative = it->second;

        full_jacobian.push_back(derivative);
        full_jacobian_sparsity.push_back(i + this->use_fortran);
    }
    return;
}
