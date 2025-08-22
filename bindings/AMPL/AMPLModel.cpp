// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <array>
#include <cassert>
#include <stdexcept>
#include "AMPLModel.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Logger.hpp"
#include "tools/Infinity.hpp"
#include "options/Options.hpp"
#include "symbolic/Concatenation.hpp"
#include "Uno.hpp"

namespace uno {
   ASL* generate_asl(const std::string &file_name) {
      ASL* asl = ASL_alloc(ASL_read_pfgh);
      FILE* nl = jac0dim_ASL(asl, file_name.data(), static_cast<int>(file_name.size()));
      // indices start at 0
      asl->i.Fortran_ = 0;

      int n_discrete = asl->i.nbv_ + asl->i.niv_ + asl->i.nlvbi_ + asl->i.nlvci_ + asl->i.nlvoi_;
      if (0 < n_discrete) {
         throw std::runtime_error("Error: " + std::to_string(n_discrete) + " variables are discrete, which Uno cannot handle");
      }

      // preallocate initial primal and dual solutions
      asl->i.X0_ = static_cast<double*>(M1zapalloc_ASL(&asl->i, sizeof(double) * static_cast<size_t>(asl->i.n_var_)));
      asl->i.pi0_ = static_cast<double*>(M1zapalloc_ASL(&asl->i, sizeof(double) * static_cast<size_t>(asl->i.n_con_)));

      // read the file_name.nl file
      pfgh_read_ASL(asl, nl, ASL_findgroups | ASL_return_read_err);

      // set the sign convention for the Lagrangian: Uno uses ∇²L(x, y) = ∇²f(x) - \sum_j y_j ∇²c_j(x)
      fint error_flag{0};
      lagscale_ASL(asl, -1., &error_flag);
      return asl;
   }

   // generate the ASL object and call the private constructor
   AMPLModel::AMPLModel(const std::string& file_name, const Options& options) : AMPLModel(file_name, generate_asl(file_name), options) {
   }

   AMPLModel::AMPLModel(const std::string& file_name, ASL* asl, const Options& options) :
         Model(file_name, static_cast<size_t>(asl->i.n_var_), static_cast<size_t>(asl->i.n_con_), (asl->i.objtype_[0] == 1) ? -1. : 1.),
         asl(asl),
         write_solution_to_file(options.get_bool("AMPL_write_solution_to_file")),
         // AMPL orders the constraints based on the function type: nonlinear first (nlc of them), then linear
         linear_constraints(static_cast<size_t>(this->asl->i.nlc_), this->number_constraints),
         equality_constraints_collection(this->equality_constraints),
         inequality_constraints_collection(this->inequality_constraints) {
      // Jacobian storage: use goff fields of struct cgrad
      this->asl->i.congrd_mode = 2;

      // detect fix variables
      this->fixed_variables.reserve(this->number_variables);
      this->partition_variables();

      // partition equality/inequality constraints
      this->equality_constraints.reserve(this->number_constraints);
      this->inequality_constraints.reserve(this->number_constraints);
      this->partition_constraints();

      // compute sparsity pattern and number of nonzeros of Lagrangian Hessian
      this->compute_lagrangian_hessian_sparsity();
   }

   AMPLModel::~AMPLModel() {
      ASL_free(&this->asl);
   }

   bool AMPLModel::has_implicit_hessian_representation() const {
      // As long as we use the ASL library ("solvers"), we need to form the explicit Hessian
      // The reason is that the ASL Hessian representation changes as soon as trial
      // iterates are evaluated. The variant "solvers2" should address the issue.
      return false;
   }

   bool AMPLModel::has_explicit_hessian_representation() const {
      return true;
   }

   double AMPLModel::evaluate_objective(const Vector<double>& x) const {
      fint error_flag = 0;
      double result = this->objective_sign * (*(this->asl)->p.Objval)(this->asl, 0, const_cast<double*>(x.data()), &error_flag);
      if (0 < error_flag) {
         throw FunctionEvaluationError();
      }
      return result;
   }

   /*
   double AMPLModel::evaluate_constraint(int j, const std::vector<double>& x) const {
      fint error_flag = 0;
      double result = (*(this->asl)->p.Conival)(this->asl_, j, const_cast<double*>(x.data()), &error_flag);
      if (0 < error_flag) {
         throw FunctionNumericalError();
      }
      return result;
   }
   */

   void AMPLModel::evaluate_constraints(const Vector<double>& x, std::vector<double>& constraints) const {
      fint error_flag = 0;
      (*(this->asl)->p.Conval)(this->asl, const_cast<double*>(x.data()), constraints.data(), &error_flag);
      if (0 < error_flag) {
         throw FunctionEvaluationError();
      }
   }

   // dense gradient
   void AMPLModel::evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const {
      fint error_flag = 0;
      (*(this->asl)->p.Objgrd)(this->asl, 0, const_cast<double*>(x.data()), gradient.data(), &error_flag);
      if (0 < error_flag) {
         throw GradientEvaluationError();
      }
      gradient.scale(this->objective_sign);
   }

   void AMPLModel::compute_constraint_jacobian_sparsity(int* row_indices, int* column_indices, int solver_indexing,
         MatrixOrder matrix_order) const {
      // by default, AMPLModel computes the Jacobian in a column-wise order (variable by variable).
      // if a row-wise evaluation is wished, modify the goff fields of the ASL "Cgrad" structures
      if (matrix_order == MatrixOrder::ROW_MAJOR) {
         int nonzero_index = 0;
         for (size_t constraint_index: Range(this->number_constraints)) {
            cgrad* constraint_gradient = this->asl->i.Cgrad_[constraint_index];
            while (constraint_gradient != nullptr) {
               constraint_gradient->goff = nonzero_index;
               ++nonzero_index;
               constraint_gradient = constraint_gradient->next;
            }
         }
      }

      for (size_t constraint_index: Range(this->number_constraints)) {
         cgrad* constraint_gradient = this->asl->i.Cgrad_[constraint_index];
         while (constraint_gradient != nullptr) {
            const int variable_index = constraint_gradient->varno;
            // at the moment, the Jacobian is stored column-wise (that is, ordered by variables)
            row_indices[constraint_gradient->goff] = static_cast<int>(constraint_index) + solver_indexing;
            column_indices[constraint_gradient->goff] = variable_index + solver_indexing;
            constraint_gradient = constraint_gradient->next;
         }
      }
   }

   void AMPLModel::compute_hessian_sparsity(int* row_indices, int* column_indices, int solver_indexing) const {
      const fint* asl_column_start = this->asl->i.sputinfo_->hcolstarts;
      const fint* asl_row_index = this->asl->i.sputinfo_->hrownos;
      size_t current_index = 0;
      for (size_t column_index: Range(this->number_variables)) {
         for (size_t nonzero_index: Range(static_cast<size_t>(asl_column_start[column_index]), static_cast<size_t>(asl_column_start[column_index + 1]))) {
            const int row_index = asl_row_index[nonzero_index];
            row_indices[current_index] = row_index + solver_indexing;
            column_indices[current_index] = static_cast<int>(column_index) + solver_indexing;
            ++current_index;
         }
      }
   }

   void AMPLModel::evaluate_constraint_jacobian(const Vector<double>& x, double* jacobian_values) const {
      fint error_flag = 0;
      (*(this->asl)->p.Jacval)(this->asl, const_cast<double*>(x.data()), jacobian_values, &error_flag);
      if (0 < error_flag) {
         throw GradientEvaluationError();
      }
   }

   // register the vector of variables
   //(*(this->asl)->p.Xknown)(this->asl, const_cast<double*>(x.data()), nullptr);
   // unregister the vector of variables
   //this->asl->i.x_known = 0;

   void AMPLModel::evaluate_lagrangian_hessian(const Vector<double>& /*x*/, double objective_multiplier, const Vector<double>& multipliers,
         double* hessian_values) const {
      objective_multiplier *= this->objective_sign;
      (*(this->asl)->p.Sphes)(this->asl, nullptr, hessian_values, -1, &objective_multiplier,
         const_cast<double*>(multipliers.data()));
   }

   void AMPLModel::compute_hessian_vector_product(const double* vector, double objective_multiplier, const Vector<double>& multipliers,
         double* result) const {
      // scale by the objective sign
      objective_multiplier *= this->objective_sign;

      // compute the Hessian-vector product
      (this->asl->p.Hvcomp)(this->asl, result, const_cast<double*>(vector), -1, &objective_multiplier,
         const_cast<double*>(multipliers.data()));
   }

   double AMPLModel::variable_lower_bound(size_t variable_index) const {
      return (this->asl->i.LUv_ != nullptr) ? this->asl->i.LUv_[2*variable_index] : -INF<double>;
   }

   double AMPLModel::variable_upper_bound(size_t variable_index) const {
      return (this->asl->i.LUv_ != nullptr) ? this->asl->i.LUv_[2*variable_index + 1] : INF<double>;
   }

   const SparseVector<size_t>& AMPLModel::get_slacks() const {
      return this->slacks;
   }

   const Vector<size_t>& AMPLModel::get_fixed_variables() const {
      return this->fixed_variables;
   }

   double AMPLModel::constraint_lower_bound(size_t constraint_index) const {
      return (this->asl->i.LUrhs_ != nullptr) ? this->asl->i.LUrhs_[2*constraint_index] : -INF<double>;
   }

   double AMPLModel::constraint_upper_bound(size_t constraint_index) const {
      return (this->asl->i.LUrhs_ != nullptr) ? this->asl->i.LUrhs_[2*constraint_index + 1] : INF<double>;
   }

   const Collection<size_t>& AMPLModel::get_equality_constraints() const {
      return this->equality_constraints_collection;
   }

   const Collection<size_t>& AMPLModel::get_inequality_constraints() const {
      return this->inequality_constraints_collection;
   }

   const Collection<size_t>& AMPLModel::get_linear_constraints() const {
      return this->linear_constraints;
   }

   // initial primal point
   void AMPLModel::initial_primal_point(Vector<double>& x) const {
      assert(x.size() >= this->number_variables);
      std::copy_n(this->asl->i.X0_, this->number_variables, x.begin());
   }

   // initial dual point
   void AMPLModel::initial_dual_point(Vector<double>& multipliers) const {
      assert(multipliers.size() >= this->number_constraints);
      std::copy_n(this->asl->i.pi0_, this->number_constraints, multipliers.begin());
   }

   void AMPLModel::postprocess_solution(Iterate& iterate, IterateStatus iterate_status) const {
      if (this->write_solution_to_file) {
         // write the primal-dual solution and status into a *.sol file
         this->asl->p.solve_code_ = 400; // limit
         if (iterate_status == IterateStatus::FEASIBLE_KKT_POINT) {
            this->asl->p.solve_code_ = 0;
         }
         if (iterate_status == IterateStatus::FEASIBLE_SMALL_STEP) {
            this->asl->p.solve_code_ = 100;
         }
         else if (iterate_status == IterateStatus::INFEASIBLE_STATIONARY_POINT) {
            this->asl->p.solve_code_ = 200;
         }
         else if (iterate_status == IterateStatus::UNBOUNDED) {
            this->asl->p.solve_code_ = 300;
         }
         else if (iterate_status == IterateStatus::INFEASIBLE_SMALL_STEP) {
            this->asl->p.solve_code_ = 500;
         }

         // flip the signs of the multipliers and the objective if we maximize
         // note: due to the different sign convention for the Lagrangian between ASL and Uno,
         // we need to flip the signs of the constraint multipliers when minimizing
         iterate.multipliers.constraints *= -this->objective_sign;
         iterate.multipliers.lower_bounds *= this->objective_sign;
         iterate.multipliers.upper_bounds *= this->objective_sign;
         iterate.evaluations.objective *= this->objective_sign;

         // include the bound duals in the .sol file, using suffixes
         SufDecl lower_bound_suffix{const_cast<char*>("lower_bound_duals"), nullptr, ASL_Sufkind_var | ASL_Sufkind_real, 0};
         SufDecl upper_bound_suffix{const_cast<char*>("upper_bound_duals"), nullptr, ASL_Sufkind_var | ASL_Sufkind_real, 0};
         std::array<SufDecl, 2> suffixes{lower_bound_suffix, upper_bound_suffix};
         suf_declare_ASL(this->asl, suffixes.data(), suffixes.size());
         suf_rput_ASL(this->asl, "lower_bound_duals", ASL_Sufkind_var, iterate.multipliers.lower_bounds.data());
         suf_rput_ASL(this->asl, "upper_bound_duals", ASL_Sufkind_var, iterate.multipliers.upper_bounds.data());

         Option_Info option_info{};
         option_info.wantsol = 9; // write the solution without printing the message to stdout
         std::string message = "Uno ";
         message.append(Uno::current_version()).append(": ").append(iterate_status_to_message(iterate_status));
         write_sol_ASL(this->asl, message.data(), iterate.primals.data(), iterate.multipliers.constraints.data(), &option_info);

         // flip back the signs of the multipliers and the objective back if we maximize
         iterate.multipliers.constraints *= -this->objective_sign;
         iterate.multipliers.lower_bounds *= this->objective_sign;
         iterate.multipliers.upper_bounds *= this->objective_sign;
         iterate.evaluations.objective *= this->objective_sign;
      }
   }

   size_t AMPLModel::number_jacobian_nonzeros() const {
      return static_cast<size_t>(this->asl->i.nzc_);
   }

   size_t AMPLModel::number_hessian_nonzeros() const {
      return this->number_asl_hessian_nonzeros;
   }

   void AMPLModel::partition_variables() {
      for (size_t variable_index: Range(this->number_variables)) {
         if (this->variable_lower_bound(variable_index) == this->variable_upper_bound(variable_index)) {
            WARNING << "Variable x" << variable_index << " has identical bounds\n";
            this->fixed_variables.emplace_back(variable_index);
         }
      }
   }

   void AMPLModel::partition_constraints() {
      for (size_t constraint_index: Range(this->number_constraints)) {
         const double lower_bound = this->constraint_lower_bound(constraint_index);
         const double upper_bound = this->constraint_upper_bound(constraint_index);
         if (lower_bound == upper_bound) {
            this->equality_constraints.emplace_back(constraint_index);
         }
         else if (!is_finite(lower_bound) && !is_finite(upper_bound)) {
            WARNING << "Constraint c" << constraint_index << " has no bounds\n";
            // count the constraint as inequality to avoid reindexing of the constraints
            this->inequality_constraints.emplace_back(constraint_index);
         }
         else {
            this->inequality_constraints.emplace_back(constraint_index);
         }
      }
   }

   void AMPLModel::compute_lagrangian_hessian_sparsity() {
      // compute the maximum number of nonzero elements, provided that all multipliers are non-zero
      // int (*Sphset) (ASL*, SputInfo**, int nobj, int ow, int y, int uptri);
      // store in lower-triangular part
      constexpr int triangular = 2;
      this->number_asl_hessian_nonzeros = static_cast<size_t>((*(this->asl)->p.Sphset)(this->asl, nullptr, -1, 1, 1, triangular));

      // sparsity pattern
      [[maybe_unused]] const fint* asl_column_start = this->asl->i.sputinfo_->hcolstarts;
      // check that the column pointers are sorted in increasing order
      assert(in_increasing_order(asl_column_start, this->number_variables + 1) &&
         "AMPLModel::evaluate_lagrangian_hessian: column starts are not ordered");
   }
} // namespace