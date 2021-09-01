#include <cassert>
#include "AMPLModel.hpp"
#include "linear_algebra/Vector.hpp"
#include "tools/Logger.hpp"

/* TODO: avoid using implicit AMPL macros */

ASL* generate_asl(std::string file_name) {
   ASL* asl = ASL_alloc(ASL_read_pfgh);
   FILE* nl = jac0dim_ASL(asl, const_cast<char*> (file_name.data()), (fint) file_name.size());
   /* indices start at 0 */
   asl->i.Fortran_ = 0;

   fint n_discrete = nlogv + niv + nlvbi + nlvci + nlvoi;
   if (0 < n_discrete) {
      WARNING << "Ignoring integrality of " << n_discrete << " variables.\n";
      asl->i.need_nl_ = 0;
   }

   /* preallocate initial primal and dual solutions */
   asl->i.X0_ = (double*) M1zapalloc(sizeof(double) * n_var);
   asl->i.pi0_ = (double*) M1zapalloc(sizeof(double) * n_con);

   /* read the file_name.nl file */
   pfgh_read_ASL(asl, nl, ASL_findgroups);

   return asl;
}

// generate the ASL object and call the private constructor
AMPLModel::AMPLModel(std::string file_name) : AMPLModel(file_name, generate_asl(file_name)) {
}

AMPLModel::AMPLModel(std::string file_name, ASL* asl) : Problem(file_name, asl->i.n_var_, asl->i.n_con_, NONLINEAR),
// asl->i.nlc_ + asl->i.nlo_ > 0),
      asl_(asl), ampl_tmp_gradient_(asl->i.n_var_) {
   this->asl_->i.congrd_mode = 0;

   /* dimensions */
   this->objective_sign = (this->asl_->i.objtype_[0] == 1) ? -1. : 1.;

   /* variables */
   this->generate_variables_();

   /* objective function */
   this->initialize_objective_();

   /* constraints */
   this->generate_constraints_();
   set_function_types_(file_name);

   /* Lagrangian Hessian */
   this->initialize_lagrangian_hessian_();
}

AMPLModel::~AMPLModel() {
   ASL_free((ASL * *) & this->asl_);
}

void AMPLModel::generate_variables_() {
   for (size_t i = 0; i < this->number_variables; i++) {
      this->variables_names[i] = var_name_ASL(this->asl_, i);
      double lb = (this->asl_->i.LUv_ != nullptr) ? this->asl_->i.LUv_[2 * i] : -INFINITY;
      double ub = (this->asl_->i.LUv_ != nullptr) ? this->asl_->i.LUv_[2 * i + 1] : INFINITY;
      if (lb == ub) {
         WARNING << "Variable x" << i << " has identical bounds\n";
      }
      this->variables_bounds[i] = {lb, ub};
   }
   this->determine_bounds_types(this->variables_bounds, this->variable_status);
}

double AMPLModel::evaluate_objective(const std::vector<double>& x) const {
   int nerror = 0;
   double result = this->objective_sign * (*(this->asl_)->p.Objval)(this->asl_, 0, (double*) x.data(), &nerror);
   if (0 < nerror) {
      throw FunctionNumericalError();
   }
   return result;
}

/* sparse gradient */
void AMPLModel::evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const {
   /* compute the AMPL gradient (always in dense format) */
   int nerror = 0;
   (*(this->asl_)->p.Objgrd)(this->asl_, 0, (double*) x.data(), (double*) this->ampl_tmp_gradient_.data(), &nerror);
   if (0 < nerror) {
      throw GradientNumericalError();
   }

   /* partial derivatives in same order as variables in this->asl_->i.Ograd_[0] */
   ograd* ampl_variables_tmp = this->asl_->i.Ograd_[0];
   while (ampl_variables_tmp != nullptr) {
      double partial_derivative = this->ampl_tmp_gradient_[ampl_variables_tmp->varno];
      /* if maximization, take the opposite */
      if (this->objective_sign < 0.) {
         partial_derivative = -partial_derivative;
      }
      gradient.insert(ampl_variables_tmp->varno, partial_derivative);
      //gradient[ampl_variables_tmp->varno] = partial_derivative;
      ampl_variables_tmp = ampl_variables_tmp->next;
   }
}

void AMPLModel::initialize_objective_() {
   this->objective_name = obj_name_ASL(this->asl_, 0);
   //this->create_objective_variables_(this->asl_->i.Ograd_[0]);
}

double AMPLModel::evaluate_constraint(int j, const std::vector<double>& x) const {
   int nerror = 0;
   double result = (*(this->asl_)->p.Conival)(this->asl_, j, (double *) x.data(), &nerror);
   if (0 < nerror) {
      throw FunctionNumericalError();
   }
   return result;
}

void AMPLModel::evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const {
   int nerror = 0;
   (*(this->asl_)->p.Conval)(this->asl_, (double *) x.data(), constraints.data(), &nerror);
   if (0 < nerror) {
      throw FunctionNumericalError();
   }
}

/* sparse gradient */
void AMPLModel::constraint_gradient(const std::vector<double>& x, int j, SparseVector<double>& gradient) const {
   const int congrd_mode_backup = this->asl_->i.congrd_mode;
   this->asl_->i.congrd_mode = 1; // sparse computation

   /* compute the AMPL gradient */
   int nerror = 0;
   (*(this->asl_)->p.Congrd)(this->asl_, j, (double*) x.data(), (double*) this->ampl_tmp_gradient_.data(), &nerror);
   if (0 < nerror) {
      throw GradientNumericalError();
   }

   /* partial derivatives in ampl_gradient in same order as variables in this->asl_->i.Cgrad_[j] */
   cgrad* ampl_variables_tmp = this->asl_->i.Cgrad_[j];
   int cpt = 0;
   while (ampl_variables_tmp != nullptr) {
      /* keep the gradient sparse */
      if (this->ampl_tmp_gradient_[cpt] != 0.) {
         gradient.insert(ampl_variables_tmp->varno, this->ampl_tmp_gradient_[cpt]);
      }
      ampl_variables_tmp = ampl_variables_tmp->next;
      cpt++;
   }
   this->asl_->i.congrd_mode = congrd_mode_backup;
}

void AMPLModel::constraints_jacobian(const std::vector<double>& x, std::vector<SparseVector<double>>& constraints_jacobian) const {
   for (size_t j = 0; j < this->number_constraints; j++) {
      this->constraint_gradient(x, j, constraints_jacobian[j]);
   }
}

void AMPLModel::generate_constraints_() {
   for (size_t j = 0; j < this->number_constraints; j++) {
      this->constraint_name[j] = con_name_ASL(this->asl_, j);
      double lb = (this->asl_->i.LUrhs_ != nullptr) ? this->asl_->i.LUrhs_[2 * j] : -INFINITY;
      double ub = (this->asl_->i.LUrhs_ != nullptr) ? this->asl_->i.LUrhs_[2 * j + 1] : INFINITY;
      this->constraint_bounds[j] = {lb, ub};
   }
   this->determine_bounds_types(this->constraint_bounds, this->constraint_status);
   this->determine_constraints();
}

void AMPLModel::set_function_types_(std::string file_name) {
   /* allocate a temporary ASL to read Hessian sparsity pattern */
   ASL* asl = ASL_alloc(ASL_read_fg);
   // char* stub = getstops(file_name, option_info);
   //if (file_name == nullptr) {
   //	usage_ASL(option_info, 1);
   //}

   FILE* nl = jac0dim(const_cast<char*> (file_name.data()), (fint) file_name.size());
   /* specific read function */
   qp_read_ASL(asl, nl, ASL_findgroups);

   fint* rowq;
   fint* colqp;
   double* delsqp;

   /* constraints */
   if ((unsigned int) asl->i.n_con_ != this->number_constraints) {
      throw std::length_error("AMPLModel.set_function_types: inconsistent number of constraints");
   }
   this->constraint_type.reserve(this->number_constraints);

   // determine the type of each constraint and objective function
   // determine if the problem is nonlinear (nonquadratic objective or nonlinear constraints)
   this->type = LINEAR;
   int current_linear_constraint = 0;
   for (size_t j = 0; j < this->number_constraints; j++) {
      fint qp = nqpcheck_ASL(asl, -(j + 1), &rowq, &colqp, &delsqp);

      if (0 < qp) {
         this->constraint_type[j] = QUADRATIC;
         this->type = NONLINEAR;
      }
      else if (qp == 0) {
         this->constraint_type[j] = LINEAR;
         this->linear_constraints[j] = current_linear_constraint;
         current_linear_constraint++;
      }
      else {
         this->constraint_type[j] = NONLINEAR;
         this->type = NONLINEAR;
      }
   }
   /* objective function */
   fint qp = nqpcheck_ASL(asl, 0, &rowq, &colqp, &delsqp);
   if (0 < qp) {
      this->objective_type = QUADRATIC;
      if (this->type == LINEAR) {
         this->type = QUADRATIC;
      }
   }
   else if (qp == 0) {
      this->objective_type = LINEAR;
   }
   else {
      this->objective_type = NONLINEAR;
      this->type = NONLINEAR;
   }
   qp_opify_ASL(asl);

   /* deallocate memory */
   ASL_free((ASL * *) & asl);
}

void AMPLModel::initialize_lagrangian_hessian_() {
   /* compute the maximum number of nonzero elements, provided that all multipliers are non-zero */
   /* fint (*Sphset) (ASL*, SputInfo**, int nobj, int ow, int y, int uptri); */
   const int objective_number = -1;
   const int upper_triangular = 1;
   this->hessian_maximum_number_nonzeros = (*(this->asl_)->p.Sphset)(this->asl_, nullptr, objective_number, 1, 1, upper_triangular);
   this->ampl_tmp_hessian.reserve(this->hessian_maximum_number_nonzeros);

   // use Lagrangian scale: in AMPL, the Lagrangian is f + lambda.g, while Uno uses f - lambda.g
   int nerror;
   lagscale_ASL(this->asl_, -1., &nerror);
}

bool are_all_zeros(const std::vector<double>& multipliers) {
   for (double xj: multipliers) {
      if (xj != 0.) {
         return false;
      }
   }
   return true;
}

size_t AMPLModel::compute_hessian_number_nonzeros(double objective_multiplier, const std::vector<double>& multipliers) const {
   /* compute the sparsity */
   const int objective_number = -1;
   const int upper_triangular = 1;
   const bool all_zeros_multipliers = are_all_zeros(multipliers);
   size_t number_non_zeros = (*(this->asl_)->p.Sphset)(this->asl_, nullptr, objective_number, (objective_multiplier > 0.),
         !all_zeros_multipliers, upper_triangular);
   return number_non_zeros;
}

void AMPLModel::lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
      CSCSymmetricMatrix& hessian) const {
   // register the vector of variables
   (*(this->asl_)->p.Xknown)(this->asl_, (double*) x.data(), 0);

   // compute the number of nonzeros
   const size_t number_non_zeros = this->fixed_hessian_sparsity ? this->hessian_maximum_number_nonzeros :
         this->compute_hessian_number_nonzeros(objective_multiplier, multipliers);
   assert(hessian.capacity >= number_non_zeros);

   // evaluate the Hessian: store the matrix in a preallocated array this->ampl_tmp_hessian
   clear(hessian.matrix);
   const int objective_number = -1;
   if (this->fixed_hessian_sparsity) {
      (*(this->asl_)->p.Sphes)(this->asl_, 0, (double*) this->ampl_tmp_hessian.data(), objective_number, &objective_multiplier,
            (double*) multipliers.data());
   }
   else {
      double* objective_multiplier_pointer = (objective_multiplier != 0.) ? &objective_multiplier : nullptr;
      bool all_zeros_multipliers = are_all_zeros(multipliers);
      (*(this->asl_)->p.Sphes)(this->asl_, 0, (double*) this->ampl_tmp_hessian.data(), objective_number, objective_multiplier_pointer,
            all_zeros_multipliers ? nullptr : (double*) multipliers.data());
   }

   // copy the nonzeros in the Hessian
   for (size_t k = 0; k < number_non_zeros; k++) {
      hessian.matrix[k] = this->ampl_tmp_hessian[k];
   }
   hessian.number_nonzeros = number_non_zeros;

   // generate the sparsity pattern in the right sparse format
   const int* ampl_column_start = this->asl_->i.sputinfo_->hcolstarts;
   const int* ampl_row_index = this->asl_->i.sputinfo_->hrownos;
   // check that the column pointers are sorted in increasing order
   assert(in_increasing_order(ampl_column_start, this->number_variables + 1) && "AMPLModel::lagrangian_hessian: column starts are not ordered");

   for (size_t k = 0; k < this->number_variables + 1; k++) {
      hessian.column_start[k] = ampl_column_start[k];
   }
   for (size_t k = 0; k < number_non_zeros; k++) {
      hessian.row_index[k] = ampl_row_index[k];
   }
   // if the subproblem has more variables (slacks, elastic, ...) than the AMPL model, rectify the sparse representation
   assert(this->number_variables <= hessian.dimension);
   for (size_t k = this->number_variables; k < hessian.dimension; k++) {
      hessian.column_start[k+1] = hessian.column_start[this->number_variables];
   }

   /* unregister the vector of variables */
   this->asl_->i.x_known = 0;
}

void AMPLModel::lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
      COOSymmetricMatrix& hessian) const {
   // register the vector of variables
   (*(this->asl_)->p.Xknown)(this->asl_, (double*) x.data(), 0);

   // compute the number of nonzeros
   const size_t number_non_zeros = this->fixed_hessian_sparsity ? this->hessian_maximum_number_nonzeros :
         this->compute_hessian_number_nonzeros(objective_multiplier, multipliers);
   assert(hessian.capacity >= number_non_zeros);

   /* evaluate the Hessian */
   clear(hessian.matrix);
   const int objective_number = -1;
   if (this->fixed_hessian_sparsity) {
      (*(this->asl_)->p.Sphes)(this->asl_, 0, (double*) hessian.matrix.data(), objective_number, &objective_multiplier,
            (double*) multipliers.data());
   }
   else {
      double* objective_multiplier_pointer = (objective_multiplier != 0.) ? &objective_multiplier : nullptr;
      bool all_zeros_multipliers = are_all_zeros(multipliers);
      (*(this->asl_)->p.Sphes)(this->asl_, 0, (double*) hessian.matrix.data(), objective_number, objective_multiplier_pointer,
            all_zeros_multipliers ? nullptr : (double*) multipliers.data());
   }
   hessian.number_nonzeros = number_non_zeros;

   // generate the sparsity pattern in the right sparse format
   const int* ampl_column_start = this->asl_->i.sputinfo_->hcolstarts;
   const int* ampl_row_index = this->asl_->i.sputinfo_->hrownos;
   int index = 0;
   for (size_t j = 0; j < this->number_variables; j++) {
      for (int k = ampl_column_start[j]; k < ampl_column_start[j + 1]; k++) {
         size_t i = ampl_row_index[k];
         hessian.row_indices[index] = i;
         hessian.column_indices[index] = j;
         index++;
      }
   }

   /* unregister the vector of variables */
   this->asl_->i.x_known = 0;
}

/* initial primal point */
void AMPLModel::set_initial_primal_point(std::vector<double>& x) {
   assert(x.size() >= this->number_variables);
   std::copy(this->asl_->i.X0_, this->asl_->i.X0_ + this->number_variables, begin(x));
}

/* initial dual point */
void AMPLModel::set_initial_dual_point(std::vector<double>& multipliers) {
   assert(multipliers.size() >= this->number_constraints);
   std::copy(this->asl_->i.pi0_, this->asl_->i.pi0_ + this->number_constraints, begin(multipliers));
}
