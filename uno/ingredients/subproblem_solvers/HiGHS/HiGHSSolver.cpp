// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "HiGHSSolver.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/Indexing.hpp"
#include "optimization/Direction.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/Logger.hpp"

namespace uno {
   HiGHSSolver::HiGHSSolver(const Options& options):
         LPSolver(), print_subproblem(options.get_bool("print_subproblem")) {
      this->highs_solver.setOptionValue("output_flag", "false");
   }

   void HiGHSSolver::initialize_memory(const Subproblem& subproblem) {
      // determine whether the subproblem has curvature. For the moment, HiGHS can only solve LPs
      if (subproblem.has_curvature()) {
         throw std::runtime_error("The subproblem has curvature. For the moment, HiGHS can only solve LPs");
      }

      this->constraints.resize(subproblem.number_constraints);
      this->linear_objective.resize(subproblem.number_variables);
      this->model.lp_.sense_ = ObjSense::kMinimize;
      this->model.lp_.offset_ = 0.;
      // the linear part of the objective is a dense vector
      this->model.lp_.col_cost_.resize(subproblem.number_variables);
      // variable bounds
      this->model.lp_.col_lower_.resize(subproblem.number_variables);
      this->model.lp_.col_upper_.resize(subproblem.number_variables);
      // constraint bounds
      this->model.lp_.row_lower_.resize(subproblem.number_constraints);
      this->model.lp_.row_upper_.resize(subproblem.number_constraints);

      // constraint Jacobian
      this->model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
      const size_t number_jacobian_nonzeros = subproblem.number_jacobian_nonzeros();
      // compute the COO sparsity pattern
      this->jacobian_row_indices.resize(number_jacobian_nonzeros);
      this->jacobian_column_indices.resize(number_jacobian_nonzeros);
      subproblem.compute_constraint_jacobian_sparsity(this->jacobian_row_indices.data(), this->jacobian_column_indices.data(),
         Indexing::C_indexing, MatrixOrder::COLUMN_MAJOR);
      // HiGHS matrix in CSC format (variable after variable)
      this->model.lp_.a_matrix_.index_.resize(number_jacobian_nonzeros); // constraint indices
      this->model.lp_.a_matrix_.start_.resize(subproblem.number_variables + 1);
      this->model.lp_.a_matrix_.value_.resize(number_jacobian_nonzeros);

      size_t current_variable = 0;
      for (size_t jacobian_nonzero_index: Range(number_jacobian_nonzeros)) {
         // constraint index
         const size_t constraint_index = this->jacobian_row_indices[jacobian_nonzero_index];
         this->model.lp_.a_matrix_.index_[jacobian_nonzero_index] = static_cast<HighsInt>(constraint_index);

         // variable index
         const size_t variable_index = this->jacobian_column_indices[jacobian_nonzero_index];
         assert(current_variable <= variable_index);
         while (current_variable < variable_index) {
            ++current_variable;
            this->model.lp_.a_matrix_.start_[current_variable] = static_cast<HighsInt>(jacobian_nonzero_index);
         }
      }
      this->model.lp_.a_matrix_.start_[subproblem.number_variables] = static_cast<HighsInt>(number_jacobian_nonzeros);

      // TODO Hessian
   }

   void HiGHSSolver::solve(Statistics& statistics, Subproblem& subproblem, const Vector<double>& /*initial_point*/,
         Direction& direction, const WarmstartInformation& warmstart_information) {
      this->set_up_subproblem(statistics, subproblem, warmstart_information);
      this->solve_subproblem(subproblem, direction);
   }

   void HiGHSSolver::compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const {
      const size_t number_constraint_jacobian_nonzeros = this->jacobian_row_indices.size();
      for (size_t nonzero_index: Range(number_constraint_jacobian_nonzeros)) {
         const size_t constraint_index = this->jacobian_row_indices[nonzero_index];
         const size_t variable_index = this->jacobian_column_indices[nonzero_index];
         const double derivative = this->model.lp_.a_matrix_.value_[nonzero_index];

         // a safeguard to make sure we take only the correct part of the Jacobian
         if (variable_index < vector.size() && constraint_index < result.size()) {
            result[constraint_index] += derivative * vector[variable_index];
         }
      }
   }

   void HiGHSSolver::compute_constraint_jacobian_transposed_vector_product(const Vector<double>& vector, Vector<double>& result) const {
      const size_t number_constraint_jacobian_nonzeros = this->jacobian_row_indices.size();
      for (size_t nonzero_index: Range(number_constraint_jacobian_nonzeros)) {
         const size_t constraint_index = this->jacobian_row_indices[nonzero_index];
         const size_t variable_index = this->jacobian_column_indices[nonzero_index];
         const double derivative = this->model.lp_.a_matrix_.value_[nonzero_index];
         assert(constraint_index < vector.size());
         assert(variable_index < result.size());

         result[variable_index] += derivative * vector[constraint_index];
      }
   }

   double HiGHSSolver::compute_hessian_quadratic_product(const Vector<double>& /*vector*/) const {
      return 0.;
   }

   void HiGHSSolver::set_up_subproblem(Statistics& /*statistics*/, const Subproblem& subproblem, const WarmstartInformation& warmstart_information) {
      this->model.lp_.num_col_ = static_cast<HighsInt>(subproblem.number_variables);
      this->model.lp_.num_row_ = static_cast<HighsInt>(subproblem.number_constraints);

      // evaluate the functions based on warmstart information
      if (warmstart_information.objective_changed) {
         subproblem.evaluate_objective_gradient(this->linear_objective);
      }
      if (warmstart_information.constraints_changed) {
         subproblem.evaluate_constraints(this->constraints);
         subproblem.evaluate_constraint_jacobian(this->model.lp_.a_matrix_.value_.data());
      }
      // evaluate the Hessian and regularize it TODO: store it in HiGHS format
      if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {
         // subproblem.compute_regularized_hessian(statistics, this->hessian);
      }

      // variable bounds
      if (warmstart_information.variable_bounds_changed) {
         subproblem.set_variables_bounds(this->model.lp_.col_lower_, this->model.lp_.col_upper_);
      }

      // constraint bounds
      if (warmstart_information.constraint_bounds_changed || warmstart_information.constraints_changed) {
         subproblem.set_constraints_bounds(this->model.lp_.row_lower_, this->model.lp_.row_upper_, this->constraints);
      }

      // copy the linear objective into the HiGHS vector
      for (size_t variable_index: Range(subproblem.number_variables)) {
         this->model.lp_.col_cost_[variable_index] = this->linear_objective[variable_index];
      }

      if (this->print_subproblem) {
         DEBUG << "LP:\n";
         DEBUG << "Linear objective part: "; print_vector(DEBUG, view(this->model.lp_.col_cost_, 0, subproblem.number_variables));
         DEBUG << "Jacobian:\n";
         DEBUG << "J = "; print_vector(DEBUG, this->model.lp_.a_matrix_.value_);
         DEBUG << "with column start: "; print_vector(DEBUG, this->model.lp_.a_matrix_.start_);
         DEBUG << "and row index: "; print_vector(DEBUG, this->model.lp_.a_matrix_.index_);
         for (size_t variable_index = 0; variable_index < subproblem.number_variables; variable_index++) {
            DEBUG << "d" << variable_index << " in [" << this->model.lp_.col_lower_[variable_index] << ", " <<
               this->model.lp_.col_upper_[variable_index] << "]\n";
         }
         for (size_t constraint_index = 0; constraint_index < subproblem.number_constraints; constraint_index++) {
            DEBUG << "linearized c" << constraint_index << " in [" << this->model.lp_.row_lower_[constraint_index] << ", " <<
               this->model.lp_.row_upper_[constraint_index]<< "]\n";
         }
      }
   }

   void HiGHSSolver::solve_subproblem(const Subproblem& subproblem, Direction& direction) {
      // solve the LP
      HighsStatus return_status = this->highs_solver.passModel(this->model);
      //assert(return_status == HighsStatus::kOk);

      return_status = this->highs_solver.run(); // solve
      DEBUG << "HiGHS status: " << static_cast<int>(return_status) << '\n';

      // if HiGHS could not optimize (e.g. because of indefinite Hessian), return an error
      if (return_status == HighsStatus::kError) {
         direction.status = SubproblemStatus::ERROR;
         return;
      }
      HighsModelStatus model_status = highs_solver.getModelStatus();
      DEBUG << "HiGHS model status: " << static_cast<int>(model_status) << '\n';

      if (model_status == HighsModelStatus::kInfeasible) {
         direction.status = SubproblemStatus::INFEASIBLE;
         return;
      }
      else if (model_status == HighsModelStatus::kUnbounded) {
         direction.status = SubproblemStatus::UNBOUNDED_PROBLEM;
         return;
      }
      
      direction.status = SubproblemStatus::OPTIMAL;
      const HighsSolution& solution = this->highs_solver.getSolution();
      // read the primal solution and bound dual solution
      for (size_t variable_index = 0; variable_index < subproblem.number_variables; variable_index++) {
         direction.primals[variable_index] = solution.col_value[variable_index];
         const double bound_multiplier = solution.col_dual[variable_index];
         if (0. < bound_multiplier) {
            direction.multipliers.lower_bounds[variable_index] = bound_multiplier;
         }
         else {
            direction.multipliers.upper_bounds[variable_index] = bound_multiplier;
         }
      }
      // read the dual solution
      for (size_t constraint_index = 0; constraint_index < subproblem.number_constraints; constraint_index++) {
         direction.multipliers.constraints[constraint_index] = solution.row_dual[constraint_index];
      }
      const HighsInfo& info = this->highs_solver.getInfo();
      direction.subproblem_objective = info.objective_function_value;
   }
} // namespace