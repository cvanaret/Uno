// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "AMPLModel.hpp"
#include "linear_algebra/Vector.hpp"
#include "tools/Logger.hpp"
#include "tools/Infinity.hpp"

ASL* generate_asl(std::string file_name) {
   ASL* asl = ASL_alloc(ASL_read_pfgh);
   FILE* nl = jac0dim_ASL(asl, file_name.data(), static_cast<int>(file_name.size()));
   // indices start at 0
   asl->i.Fortran_ = 0;

   int n_discrete = asl->i.nbv_ + asl->i.niv_ + asl->i.nlvbi_ + asl->i.nlvci_ + asl->i.nlvoi_;
   if (0 < n_discrete) {
      WARNING << "Ignoring integrality of " << n_discrete << " variables.\n";
      asl->i.need_nl_ = 0;
   }

   // preallocate initial primal and dual solutions
   asl->i.X0_ = static_cast<double*>(M1zapalloc_ASL(&asl->i, sizeof(double) * static_cast<size_t>(n_var)));
   asl->i.pi0_ = static_cast<double*>(M1zapalloc_ASL(&asl->i, sizeof(double) * static_cast<size_t>(n_con)));

   // read the file_name.nl file
   pfgh_read_ASL(asl, nl, ASL_findgroups);
   return asl;
}

// generate the ASL object and call the private constructor
AMPLModel::AMPLModel(const std::string& file_name) : AMPLModel(file_name, generate_asl(file_name)) {
}

AMPLModel::AMPLModel(const std::string& file_name, ASL* asl) :
      Model(file_name, static_cast<size_t>(asl->i.n_var_), static_cast<size_t>(asl->i.n_con_), NONLINEAR),
      asl(asl),
      // allocate vectors
      ampl_tmp_gradient(this->number_variables),
      variables_bounds(this->number_variables),
      constraint_bounds(this->number_constraints),
      variable_status(this->number_variables),
      constraint_type(this->number_constraints),
      constraint_status(this->number_constraints) {
   this->asl->i.congrd_mode = 0;

   // dimensions
   this->objective_sign = (this->asl->i.objtype_[0] == 1) ? -1. : 1.;

   // variables
   this->generate_variables();

   // constraints
   this->equality_constraints.reserve(this->number_constraints);
   this->inequality_constraints.reserve(this->number_constraints);
   this->linear_constraints.reserve(this->number_constraints);
   this->generate_constraints();
   this->set_function_types(file_name);

   // compute number of nonzeros
   this->number_objective_gradient_nonzeros = static_cast<size_t>(this->asl->i.nzo_);
   this->number_jacobian_nonzeros = static_cast<size_t>(this->asl->i.nzc_);
   this->set_number_hessian_nonzeros();
}

AMPLModel::~AMPLModel() {
   ASL_free(&this->asl);
}

void AMPLModel::generate_variables() {
   for (size_t i: Range(this->number_variables)) {
      double lb = (this->asl->i.LUv_ != nullptr) ? this->asl->i.LUv_[2 * i] : -INF<double>;
      double ub = (this->asl->i.LUv_ != nullptr) ? this->asl->i.LUv_[2 * i + 1] : INF<double>;
      if (lb == ub) {
         WARNING << "Variable x" << i << " has identical bounds\n";
      }
      this->variables_bounds[i] = {lb, ub};
   }
   Model::determine_bounds_types(this->variables_bounds, this->variable_status);
   // figure out the bounded variables
   for (size_t i: Range(this->number_variables)) {
      const BoundType status = this->get_variable_bound_type(i);
      if (status == BOUNDED_LOWER || status == BOUNDED_BOTH_SIDES) {
         this->lower_bounded_variables.push_back(i);
         if (status == BOUNDED_LOWER) {
            this->single_lower_bounded_variables.push_back(i);
         }
      }
      if (status == BOUNDED_UPPER || status == BOUNDED_BOTH_SIDES) {
         this->upper_bounded_variables.push_back(i);
         if (status == BOUNDED_UPPER) {
            this->single_upper_bounded_variables.push_back(i);
         }
      }
   }
}

double AMPLModel::evaluate_objective(const std::vector<double>& x) const {
   int nerror = 0;
   double result = this->objective_sign * (*(this->asl)->p.Objval)(this->asl, 0, const_cast<double*>(x.data()), &nerror);
   if (0 < nerror) {
      throw FunctionEvaluationError();
   }
   return result;
}

// sparse gradient
void AMPLModel::evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const {
   // compute the AMPL gradient (always in dense format)
   int nerror = 0;
   (*(this->asl)->p.Objgrd)(this->asl, 0, const_cast<double*>(x.data()), const_cast<double*>(this->ampl_tmp_gradient.data()), &nerror);
   if (0 < nerror) {
      throw GradientEvaluationError();
   }

   // partial derivatives in same order as variables in this->asl_->i.Ograd_[0]
   ograd* ampl_variables_tmp = this->asl->i.Ograd_[0];
   while (ampl_variables_tmp != nullptr) {
      const size_t index = static_cast<size_t>(ampl_variables_tmp->varno);
      // scale by the objective sign
      const double partial_derivative = this->objective_sign*this->ampl_tmp_gradient[index];
      gradient.insert(index, partial_derivative);
      ampl_variables_tmp = ampl_variables_tmp->next;
   }
}

/*
double AMPLModel::evaluate_constraint(int j, const std::vector<double>& x) const {
   int nerror = 0;
   double result = (*(this->asl)->p.Conival)(this->asl_, j, const_cast<double*>(x.data()), &nerror);
   if (0 < nerror) {
      throw FunctionNumericalError();
   }
   return result;
}
*/

void AMPLModel::evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const {
   int nerror = 0;
   (*(this->asl)->p.Conval)(this->asl, const_cast<double*>(x.data()), constraints.data(), &nerror);
   if (0 < nerror) {
      throw FunctionEvaluationError();
   }
}

// sparse gradient
void AMPLModel::evaluate_constraint_gradient(const std::vector<double>& x, size_t j, SparseVector<double>& gradient) const {
   const int congrd_mode_backup = this->asl->i.congrd_mode;
   this->asl->i.congrd_mode = 1; // sparse computation

   // compute the AMPL gradient
   int nerror = 0;
   (*(this->asl)->p.Congrd)(this->asl, static_cast<int>(j), const_cast<double*>(x.data()), const_cast<double*>(this->ampl_tmp_gradient.data()),
         &nerror);
   if (0 < nerror) {
      throw GradientEvaluationError();
   }

   // partial derivatives in ampl_gradient in same order as variables in this->asl_->i.Cgrad_[j]
   gradient.clear();
   cgrad* ampl_variables_tmp = this->asl->i.Cgrad_[j];
   size_t index = 0;
   while (ampl_variables_tmp != nullptr) {
      gradient.insert(static_cast<size_t>(ampl_variables_tmp->varno), this->ampl_tmp_gradient[index]);
      ampl_variables_tmp = ampl_variables_tmp->next;
      index++;
   }
   this->asl->i.congrd_mode = congrd_mode_backup;
}

void AMPLModel::evaluate_constraint_jacobian(const std::vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const {
   for (size_t j: Range(this->number_constraints)) {
      constraint_jacobian[j].clear();
      this->evaluate_constraint_gradient(x, j, constraint_jacobian[j]);
   }
}

void AMPLModel::set_number_hessian_nonzeros() {
   // compute the maximum number of nonzero elements, provided that all multipliers are non-zero
   // int (*Sphset) (ASL*, SputInfo**, int nobj, int ow, int y, int uptri);
   const int objective_number = -1;
   const int upper_triangular = 1;
   this->number_hessian_nonzeros = static_cast<size_t>((*(this->asl)->p.Sphset)(this->asl, nullptr, objective_number, 1, 1, upper_triangular));
   this->ampl_tmp_hessian.reserve(this->number_hessian_nonzeros);

   // use Lagrangian scale: in AMPL, the Lagrangian is f + lambda.g, while Uno uses f - lambda.g
   int nerror{};
   lagscale_ASL(this->asl, -1., &nerror);
}

size_t AMPLModel::get_number_objective_gradient_nonzeros() const {
   return this->number_objective_gradient_nonzeros;
}

size_t AMPLModel::get_number_jacobian_nonzeros() const {
   return this->number_jacobian_nonzeros;
}

size_t AMPLModel::get_number_hessian_nonzeros() const {
   return this->number_hessian_nonzeros;
}

bool are_all_zeros(const std::vector<double>& multipliers) {
   return std::all_of(multipliers.cbegin(), multipliers.cend(), [](double xj) {
      return xj == 0.;
   });
}

size_t AMPLModel::compute_hessian_number_nonzeros(double objective_multiplier, const std::vector<double>& multipliers) const {
   // compute the sparsity
   const int objective_number = -1;
   const int upper_triangular = 1;
   const bool all_zeros_multipliers = are_all_zeros(multipliers);
   int number_nonzeros = (*(this->asl)->p.Sphset)(this->asl, nullptr, objective_number, (objective_multiplier != 0.),
         not all_zeros_multipliers, upper_triangular);
   return static_cast<size_t>(number_nonzeros);
}

void AMPLModel::evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
      SymmetricMatrix<double>& hessian) const {
   // register the vector of variables
   (*(this->asl)->p.Xknown)(this->asl, const_cast<double*>(x.data()), nullptr);

   // scale by the objective sign
   objective_multiplier *= this->objective_sign;

   // compute the number of nonzeros
   [[maybe_unused]] const size_t number_nonzeros = this->fixed_hessian_sparsity ? this->number_hessian_nonzeros :
                                                   this->compute_hessian_number_nonzeros(objective_multiplier, multipliers);
   assert(hessian.capacity >= number_nonzeros);

   // evaluate the Hessian: store the matrix in a preallocated array this->ampl_tmp_hessian
   const int objective_number = -1;
   if (this->fixed_hessian_sparsity) {
      (*(this->asl)->p.Sphes)(this->asl, nullptr, const_cast<double*>(this->ampl_tmp_hessian.data()), objective_number, &objective_multiplier,
            const_cast<double*>(multipliers.data()));
   }
   else {
      double* objective_multiplier_pointer = (objective_multiplier != 0.) ? &objective_multiplier : nullptr;
      bool all_zeros_multipliers = are_all_zeros(multipliers);
      (*(this->asl)->p.Sphes)(this->asl, nullptr, const_cast<double*>(this->ampl_tmp_hessian.data()), objective_number, objective_multiplier_pointer,
            all_zeros_multipliers ? nullptr : const_cast<double*>(multipliers.data()));
   }

   // generate the sparsity pattern in the right sparse format
   const int* ampl_column_start = this->asl->i.sputinfo_->hcolstarts;
   const int* ampl_row_index = this->asl->i.sputinfo_->hrownos;
   // check that the column pointers are sorted in increasing order
   assert(in_increasing_order(ampl_column_start, this->number_variables + 1) && "AMPLModel::evaluate_lagrangian_hessian: column starts are not ordered");

   // copy the nonzeros in the Hessian
   hessian.reset();
   for (size_t j: Range(this->number_variables)) {
      for (size_t k: Range(static_cast<size_t>(ampl_column_start[j]), static_cast<size_t>(ampl_column_start[j + 1]))) {
         const size_t i = static_cast<size_t>(ampl_row_index[k]);
         const double entry = this->ampl_tmp_hessian[k];
         hessian.insert(entry, i, j);
      }
      hessian.finalize_column(j);
   }
   // unregister the vector of variables
   this->asl->i.x_known = 0;
}

double AMPLModel::get_variable_lower_bound(size_t i) const {
   return this->variables_bounds[i].lb;
}

double AMPLModel::get_variable_upper_bound(size_t i) const {
   return this->variables_bounds[i].ub;
}

BoundType AMPLModel::get_variable_bound_type(size_t i) const {
   return this->variable_status[i];
}

double AMPLModel::get_constraint_lower_bound(size_t j) const {
   return this->constraint_bounds[j].lb;
}

double AMPLModel::get_constraint_upper_bound(size_t j) const {
   return this->constraint_bounds[j].ub;
}

FunctionType AMPLModel::get_constraint_type(size_t j) const {
   return this->constraint_type[j];
}

BoundType AMPLModel::get_constraint_bound_type(size_t j) const {
   return this->constraint_status[j];
}

// initial primal point
void AMPLModel::get_initial_primal_point(std::vector<double>& x) const {
   assert(x.size() >= this->number_variables);
   std::copy(this->asl->i.X0_, this->asl->i.X0_ + this->number_variables, begin(x));
}

// initial dual point
void AMPLModel::get_initial_dual_point(std::vector<double>& multipliers) const {
   assert(multipliers.size() >= this->number_constraints);
   std::copy(this->asl->i.pi0_, this->asl->i.pi0_ + this->number_constraints, begin(multipliers));
}

void AMPLModel::postprocess_solution(Iterate& /*iterate*/, TerminationStatus /*termination_status*/) const {
   // do nothing
}

const std::vector<size_t>& AMPLModel::get_linear_constraints() const {
   return this->linear_constraints;
}

void AMPLModel::generate_constraints() {
   for (size_t j: Range(this->number_constraints)) {
      double lb = (this->asl->i.LUrhs_ != nullptr) ? this->asl->i.LUrhs_[2 * j] : -INF<double>;
      double ub = (this->asl->i.LUrhs_ != nullptr) ? this->asl->i.LUrhs_[2 * j + 1] : INF<double>;
      this->constraint_bounds[j] = {lb, ub};
   }
   Model::determine_bounds_types(this->constraint_bounds, this->constraint_status);
   this->determine_constraints();
}

void AMPLModel::set_function_types(std::string file_name) {
   // allocate a temporary ASL to read Hessian sparsity pattern
   ASL* asl_fg = ASL_alloc(ASL_read_fg);
   // char* stub = getstops(file_name, option_info);
   //if (file_name == nullptr) {
   //	usage_ASL(option_info, 1);
   //}

   FILE* nl = jac0dim_ASL(asl_fg, file_name.data(), static_cast<int>(file_name.size()));
   // specific read function
   qp_read_ASL(asl_fg, nl, ASL_findgroups);

   // constraints
   if (asl_fg->i.n_con_ != static_cast<int>(this->number_constraints)) {
      throw std::length_error("AMPLModel.set_function_types: inconsistent number of constraints");
   }
   this->constraint_type.reserve(this->number_constraints);

   // determine the type of each constraint and objective function
   // determine if the problem is nonlinear (non-quadratic objective or nonlinear constraints)
   this->problem_type = LINEAR;
   int* rowq;
   int* colqp;
   double* delsqp;
   for (size_t j: Range(this->number_constraints)) {
      int qp = nqpcheck_ASL(asl_fg, static_cast<int>(-(j + 1)), &rowq, &colqp, &delsqp);

      if (0 < qp) {
         this->constraint_type[j] = QUADRATIC;
         this->problem_type = NONLINEAR;
      }
      else if (qp == 0) {
         this->constraint_type[j] = LINEAR;
         this->linear_constraints.push_back(j);
      }
      else {
         this->constraint_type[j] = NONLINEAR;
         this->problem_type = NONLINEAR;
      }
   }
   // objective function
   int qp = nqpcheck_ASL(asl_fg, 0, &rowq, &colqp, &delsqp);
   if (0 < qp) {
      if (this->problem_type == LINEAR) {
         this->problem_type = QUADRATIC;
      }
   }
   else if (qp != 0) {
      this->problem_type = NONLINEAR;
   }
   qp_opify_ASL(asl_fg);

   // deallocate memory
   ASL_free(&asl_fg);
}
