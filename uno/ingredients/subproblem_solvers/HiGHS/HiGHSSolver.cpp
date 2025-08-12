// Copyright (c) 2024-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "HiGHSSolver.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/Indexing.hpp"
#include "optimization/Direction.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"

namespace uno {
   HiGHSSolver::HiGHSSolver(const Options& options):
         QPSolver(), print_subproblem(options.get_bool("print_subproblem")) {
      this->highs_solver.setOptionValue("output_flag", "false");
   }

   void HiGHSSolver::initialize_memory(const Subproblem& subproblem) {
      this->evaluation_space.model.lp_.num_col_ = static_cast<HighsInt>(subproblem.number_variables);
      this->evaluation_space.model.lp_.num_row_ = static_cast<HighsInt>(subproblem.number_constraints);

      // determine whether the subproblem has curvature. For the moment, HiGHS can only solve LPs
      this->evaluation_space.constraints.resize(subproblem.number_constraints);
      this->evaluation_space.linear_objective.resize(subproblem.number_variables);
      this->evaluation_space.model.lp_.sense_ = ObjSense::kMinimize;
      this->evaluation_space.model.lp_.offset_ = 0.;
      // the linear part of the objective is a dense vector
      this->evaluation_space.model.lp_.col_cost_.resize(subproblem.number_variables);
      // variable bounds
      this->evaluation_space.model.lp_.col_lower_.resize(subproblem.number_variables);
      this->evaluation_space.model.lp_.col_upper_.resize(subproblem.number_variables);
      // constraint bounds
      this->evaluation_space.model.lp_.row_lower_.resize(subproblem.number_constraints);
      this->evaluation_space.model.lp_.row_upper_.resize(subproblem.number_constraints);

      // column-wise constraint Jacobian
      this->evaluation_space.model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
      const size_t number_jacobian_nonzeros = subproblem.number_jacobian_nonzeros();
      // compute the COO sparsity pattern
      this->evaluation_space.jacobian_row_indices.resize(number_jacobian_nonzeros);
      this->evaluation_space.jacobian_column_indices.resize(number_jacobian_nonzeros);
      subproblem.compute_constraint_jacobian_sparsity(this->evaluation_space.jacobian_row_indices.data(),
         this->evaluation_space.jacobian_column_indices.data(), Indexing::C_indexing, MatrixOrder::COLUMN_MAJOR);
      // HiGHS matrix in CSC format (variable after variable)
      this->evaluation_space.model.lp_.a_matrix_.index_.resize(number_jacobian_nonzeros); // constraint indices
      this->evaluation_space.model.lp_.a_matrix_.start_.resize(subproblem.number_variables + 1);
      this->evaluation_space.model.lp_.a_matrix_.value_.resize(number_jacobian_nonzeros);
      size_t current_variable = 0;
      for (size_t jacobian_nonzero_index: Range(number_jacobian_nonzeros)) {
         // constraint index is used as is
         const size_t constraint_index = this->evaluation_space.jacobian_row_indices[jacobian_nonzero_index];
         this->evaluation_space.model.lp_.a_matrix_.index_[jacobian_nonzero_index] = static_cast<HighsInt>(constraint_index);

         // variable index is used to build the pointers to the column starts
         const size_t variable_index = this->evaluation_space.jacobian_column_indices[jacobian_nonzero_index];
         assert(current_variable <= variable_index);
         while (current_variable < variable_index) {
            ++current_variable;
            this->evaluation_space.model.lp_.a_matrix_.start_[current_variable] = static_cast<HighsInt>(jacobian_nonzero_index);
         }
      }
      this->evaluation_space.model.lp_.a_matrix_.start_[subproblem.number_variables] = static_cast<HighsInt>(number_jacobian_nonzeros);

      // Lagrangian Hessian
      this->compute_hessian_sparsity(subproblem);
   }

   void HiGHSSolver::solve(Statistics& statistics, Subproblem& subproblem, const Vector<double>& /*initial_point*/,
         Direction& direction, const WarmstartInformation& warmstart_information) {
      this->set_up_subproblem(statistics, subproblem, warmstart_information);
      this->solve_subproblem(subproblem, direction);
   }

   EvaluationSpace& HiGHSSolver::get_evaluation_space() {
      return this->evaluation_space;
   }

   // protected member functions

   // column-wise lower triangular Lagrangian Hessian
   void HiGHSSolver::compute_hessian_sparsity(const Subproblem& subproblem) {
      const size_t number_regularized_hessian_nonzeros = subproblem.number_regularized_hessian_nonzeros();
      this->evaluation_space.model.hessian_.dim_ = static_cast<HighsInt>(subproblem.number_variables);
      this->evaluation_space.model.hessian_.format_ = HessianFormat::kTriangular;
      this->evaluation_space.model.hessian_.index_.resize(number_regularized_hessian_nonzeros);
      this->evaluation_space.model.hessian_.start_.resize(subproblem.number_variables + 1);
      this->evaluation_space.model.hessian_.value_.resize(number_regularized_hessian_nonzeros);

      // get the Jacobian sparsity in COO format
      this->evaluation_space.hessian_row_indices.resize(number_regularized_hessian_nonzeros);
      this->evaluation_space.hessian_column_indices.resize(number_regularized_hessian_nonzeros);
      subproblem.compute_regularized_hessian_sparsity(this->evaluation_space.hessian_row_indices.data(),
         this->evaluation_space.hessian_column_indices.data(), Indexing::C_indexing);

      // HiGHS requires a lower-triangular CSC Hessian: the entries should be in increasing column indices.
      // Since the COO format does not require this, we need to convert from COO to CSC by permutating the entries. To
      // this end, we compute a permutation vector once and for all that contains the correct ordering of terms.
      // The permutation vector is initially filled with [0, 1, ..., number_regularized_hessian_nonzeros-1]
      this->evaluation_space.permutation_vector.resize(number_regularized_hessian_nonzeros);
      std::iota(this->evaluation_space.permutation_vector.begin(), this->evaluation_space.permutation_vector.end(), 0);
      // sort the permutation vector such that the column indices are in increasing order
      // see https://stackoverflow.com/questions/17554242/how-to-obtain-the-index-permutation-after-the-sorting
      std::sort(this->evaluation_space.permutation_vector.begin(), this->evaluation_space.permutation_vector.end(),
          [&](const size_t& i, const size_t& j) {
             if (this->evaluation_space.hessian_column_indices[i] < this->evaluation_space.hessian_column_indices[j]) {
                return true;
             }
             // within a given column, have the row indices in increasing order
             else if (this->evaluation_space.hessian_column_indices[i] == this->evaluation_space.hessian_column_indices[j]) {
               return (this->evaluation_space.hessian_row_indices[i] < this->evaluation_space.hessian_row_indices[j]);
             }
             return false;
          }
      );

      // copy the COO format into HiGHS' CSC format
      this->evaluation_space.model.hessian_.start_[0] = 0;
      size_t current_column = 0;
      for (size_t hessian_nonzero_index: Range(number_regularized_hessian_nonzeros)) {
         const size_t permutated_nonzero_index = this->evaluation_space.permutation_vector[hessian_nonzero_index];
         // row index
         const size_t row_index = this->evaluation_space.hessian_row_indices[permutated_nonzero_index];
         this->evaluation_space.model.hessian_.index_[hessian_nonzero_index] = static_cast<HighsInt>(row_index);

         // column index
         const size_t column_index = this->evaluation_space.hessian_column_indices[permutated_nonzero_index];
         assert(current_column <= column_index);
         while (current_column < column_index) {
            ++current_column;
            this->evaluation_space.model.hessian_.start_[current_column] = static_cast<HighsInt>(hessian_nonzero_index);
         }
      }
      // fill the remaining empty columns
      while (current_column < subproblem.number_variables) {
         ++current_column;
         this->evaluation_space.model.hessian_.start_[current_column] = static_cast<HighsInt>(number_regularized_hessian_nonzeros);
      }

      // the Hessian will be evaluated in this vector, and copied with permutation into this->evaluation_space.model.hessian_.value_
      this->evaluation_space.hessian_values.resize(number_regularized_hessian_nonzeros);
   }

   void HiGHSSolver::set_up_subproblem(Statistics& statistics, const Subproblem& subproblem, const WarmstartInformation& warmstart_information) {
      this->evaluation_space.evaluate_functions(statistics, subproblem, warmstart_information);

      // variable bounds
      if (warmstart_information.variable_bounds_changed) {
         subproblem.set_variables_bounds(this->evaluation_space.model.lp_.col_lower_, this->evaluation_space.model.lp_.col_upper_);
      }

      // constraint bounds
      if (warmstart_information.constraint_bounds_changed || warmstart_information.constraints_changed) {
         subproblem.set_constraints_bounds(this->evaluation_space.model.lp_.row_lower_, this->evaluation_space.model.lp_.row_upper_, this->evaluation_space.constraints);
      }

      if (this->print_subproblem) {
         DEBUG << "Subproblem:\n";
         DEBUG << "Hessian: "; print_vector(DEBUG, this->evaluation_space.model.hessian_.value_);
         DEBUG << "Linear objective part: "; print_vector(DEBUG, this->evaluation_space.model.lp_.col_cost_);
         DEBUG << "Jacobian: "; print_vector(DEBUG, this->evaluation_space.model.lp_.a_matrix_.value_);
         // DEBUG << "with column start: "; print_vector(DEBUG, this->evaluation_space.model.lp_.a_matrix_.start_);
         // DEBUG << "and row index: "; print_vector(DEBUG, this->evaluation_space.model.lp_.a_matrix_.index_);
         for (size_t variable_index = 0; variable_index < subproblem.number_variables; variable_index++) {
            DEBUG << "d" << variable_index << " in [" << this->evaluation_space.model.lp_.col_lower_[variable_index] << ", " <<
               this->evaluation_space.model.lp_.col_upper_[variable_index] << "]\n";
         }
         for (size_t constraint_index = 0; constraint_index < subproblem.number_constraints; constraint_index++) {
            DEBUG << "linearized c" << constraint_index << " in [" << this->evaluation_space.model.lp_.row_lower_[constraint_index] << ", " <<
               this->evaluation_space.model.lp_.row_upper_[constraint_index]<< "]\n";
         }
      }
   }

   void HiGHSSolver::solve_subproblem(const Subproblem& subproblem, Direction& direction) {
      // solve the subproblem
      HighsStatus return_status = this->highs_solver.passModel(this->evaluation_space.model);
      //assert(return_status == HighsStatus::kOk);

      DEBUG2 << "Running HiGHS\n";
      return_status = this->highs_solver.run(); // solve
      DEBUG2 << "Ran HiGHS\n";
      DEBUG << "HiGHS status: " << static_cast<int>(return_status) << '\n';

      // if HiGHS could not optimize (e.g. because of indefinite Hessian), return an error
      if (return_status == HighsStatus::kError) {
         throw std::runtime_error("HiGHS encountered negative curvature, which it cannot handle. Terminating.");
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