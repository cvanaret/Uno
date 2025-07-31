// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "HiGHSSolver.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/regularization_strategies/RegularizationStrategy.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/Logger.hpp"

namespace uno {
   HiGHSSolver::HiGHSSolver(const Options& options):
         LPSolver(), print_subproblem(options.get_bool("print_subproblem")) {
      this->highs_solver.setOptionValue("output_flag", "false");
   }

   void HiGHSSolver::initialize_memory(const OptimizationProblem& problem, const HessianModel& hessian_model,
         const RegularizationStrategy<double>& regularization_strategy) {
      // determine whether the subproblem has curvature. For the moment, HiGHS can only solve LPs
      const size_t number_hessian_nonzeros = problem.number_hessian_nonzeros(hessian_model);
      const bool regularize = !hessian_model.is_positive_definite() && regularization_strategy.performs_primal_regularization();
      const size_t regularization_size = problem.get_number_original_variables();
      const size_t number_regularized_hessian_nonzeros = number_hessian_nonzeros +
         (regularize ? regularization_size : 0);
      if (0 < number_regularized_hessian_nonzeros) {
         throw std::runtime_error("The subproblem has curvature. For the moment, HiGHS can only solve LPs");
      }

      this->constraints.resize(problem.number_constraints);
      this->linear_objective.resize(problem.number_variables);
      this->constraint_jacobian.resize(problem.number_constraints, problem.number_variables);
      this->model.lp_.sense_ = ObjSense::kMinimize;
      this->model.lp_.offset_ = 0.;
      // the linear part of the objective is a dense vector
      this->model.lp_.col_cost_.resize(problem.number_variables);
      // variable bounds
      this->model.lp_.col_lower_.resize(problem.number_variables);
      this->model.lp_.col_upper_.resize(problem.number_variables);
      // constraint bounds
      this->model.lp_.row_lower_.resize(problem.number_constraints);
      this->model.lp_.row_upper_.resize(problem.number_constraints);
      // constraint matrix in CSR format (each row is a constraint gradient)
      this->model.lp_.a_matrix_.format_ = MatrixFormat::kRowwise;
      this->model.lp_.a_matrix_.value_.reserve(problem.number_jacobian_nonzeros());
      this->model.lp_.a_matrix_.index_.reserve(problem.number_jacobian_nonzeros());
      this->model.lp_.a_matrix_.start_.reserve(problem.number_variables + 1);
      // TODO Hessian
      // this->hessian = SymmetricMatrix<size_t, double>("CSC", problem.number_variables, 0, 0); // TODO
   }

   void HiGHSSolver::solve(Statistics& statistics, Subproblem& subproblem, const Vector<double>& /*initial_point*/,
         Direction& direction, const WarmstartInformation& warmstart_information) {
      this->set_up_subproblem(statistics, subproblem, warmstart_information);
      this->solve_subproblem(subproblem, direction);
   }

   double HiGHSSolver::hessian_quadratic_product(const Vector<double>& /*vector*/) const {
      return 0.;
   }

   void HiGHSSolver::set_up_subproblem(Statistics& statistics, const Subproblem& subproblem, const WarmstartInformation& warmstart_information) {
      this->model.lp_.num_col_ = static_cast<HighsInt>(subproblem.number_variables);
      this->model.lp_.num_row_ = static_cast<HighsInt>(subproblem.number_constraints);

      // evaluate the functions based on warmstart information
      if (warmstart_information.objective_changed) {
         subproblem.evaluate_objective_gradient(this->linear_objective);
      }
      if (warmstart_information.constraints_changed) {
         subproblem.evaluate_constraints(this->constraints);
         subproblem.evaluate_jacobian(this->constraint_jacobian);
      }
      // evaluate the Hessian and regularize it TODO: store it in HiGHS format
      if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {
         this->hessian.reset();
         subproblem.compute_regularized_hessian(statistics, this->hessian);
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

      // copy the constraint matrix into the HiGHS structure
      this->model.lp_.a_matrix_.value_.clear();
      this->model.lp_.a_matrix_.index_.clear();
      this->model.lp_.a_matrix_.start_.clear();
      size_t number_nonzeros = 0;
      this->model.lp_.a_matrix_.start_.emplace_back(number_nonzeros);
      for (size_t constraint_index: Range(subproblem.number_constraints)) {
         for (const auto [variable_index, value]: this->constraint_jacobian[constraint_index]) {
            this->model.lp_.a_matrix_.value_.emplace_back(value);
            this->model.lp_.a_matrix_.index_.emplace_back(variable_index);
            number_nonzeros++;
         }
         this->model.lp_.a_matrix_.start_.emplace_back(number_nonzeros);
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