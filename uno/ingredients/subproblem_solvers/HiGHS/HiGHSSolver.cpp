#include <cassert>
#include "HiGHSSolver.hpp"
#include "ingredients/constraint_relaxation_strategies/OptimizationProblem.hpp"
#include "layers/SubproblemLayer.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "symbolic/VectorView.hpp"

namespace uno {
   HiGHSSolver::HiGHSSolver(const Options& options):
         LPSolver(),
         print_subproblem(options.get_bool("print_subproblem")) {
      this->highs_solver.setOptionValue("output_flag", "false");
   }

   void HiGHSSolver::initialize_memory(const OptimizationProblem& problem, const HessianModel& /*hessian_model*/,
         RegularizationStrategy<double>& /*regularization_strategy*/) {
      this->constraints.resize(problem.number_constraints);
      this->linear_objective.reserve(problem.number_objective_gradient_nonzeros());
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
      // const size_t number_hessian_nonzeros = hessian_model.number_nonzeros(problem);
   }

   void HiGHSSolver::solve(Statistics& /*statistics*/, const OptimizationProblem& problem, Iterate& current_iterate,
         const Multipliers& /*current_multipliers*/, const Vector<double>& /*initial_point*/, Direction& direction,
         SubproblemLayer& subproblem_layer, double trust_region_radius, const WarmstartInformation& warmstart_information) {
      // check if HiGHS can solve the subproblem
      const size_t number_regularized_hessian_nonzeros = problem.number_hessian_nonzeros(*subproblem_layer.hessian_model) +
         (subproblem_layer.regularization_strategy->performs_primal_regularization() ? problem.number_variables : 0);
      if (0 < number_regularized_hessian_nonzeros) {
         throw std::runtime_error("The subproblem has curvature. For the moment, HiGHS can only solve LPs");
      }

      if (this->print_subproblem) {
         DEBUG << "LP:\n";
      }
      this->set_up_subproblem(problem, current_iterate, trust_region_radius, warmstart_information);
      this->solve_subproblem(problem, direction);
   }

   double HiGHSSolver::hessian_quadratic_product(const Vector<double>& /*vector*/) const {
      return 0.;
   }

   void HiGHSSolver::set_up_subproblem(const OptimizationProblem& problem, Iterate& current_iterate, double trust_region_radius,
         const WarmstartInformation& warmstart_information) {
      this->model.lp_.num_col_ = static_cast<HighsInt>(problem.number_variables);
      this->model.lp_.num_row_ = static_cast<HighsInt>(problem.number_constraints);

      // function evaluations
      if (warmstart_information.objective_changed) {
         problem.evaluate_objective_gradient(current_iterate, this->linear_objective);
         // reset the linear part of the objective
         for (size_t variable_index: Range(problem.number_variables)) {
            this->model.lp_.col_cost_[variable_index] = 0.;
         }
         // copy the sparse vector into the dense vector
         for (const auto [variable_index, derivative]: this->linear_objective) {
            this->model.lp_.col_cost_[variable_index] = derivative;
         }
      }
      if (warmstart_information.constraints_changed) {
         problem.evaluate_constraints(current_iterate, this->constraints);
         problem.evaluate_constraint_jacobian(current_iterate, this->constraint_jacobian);
      }

      // linear part of the objective
      for (const auto [variable_index, value]: this->linear_objective) {
         this->model.lp_.col_cost_[variable_index] = value;
      }

      // constraint matrix
      this->model.lp_.a_matrix_.value_.clear();
      this->model.lp_.a_matrix_.index_.clear();
      this->model.lp_.a_matrix_.start_.clear();

      size_t number_nonzeros = 0;
      this->model.lp_.a_matrix_.start_.emplace_back(number_nonzeros);
      for (size_t constraint_index: Range(problem.number_constraints)) {
         for (const auto [variable_index, value]: this->constraint_jacobian[constraint_index]) {
            this->model.lp_.a_matrix_.value_.emplace_back(value);
            this->model.lp_.a_matrix_.index_.emplace_back(variable_index);
            number_nonzeros++;
         }
         this->model.lp_.a_matrix_.start_.emplace_back(number_nonzeros);
      }

      // variable bounds
      if (warmstart_information.variable_bounds_changed) {
         // bounds of original variables intersected with trust region
         for (size_t variable_index: Range(problem.get_number_original_variables())) {
            this->model.lp_.col_lower_[variable_index] = std::max(-trust_region_radius,
                  problem.variable_lower_bound(variable_index) - current_iterate.primals[variable_index]);
            this->model.lp_.col_upper_[variable_index] = std::min(trust_region_radius,
                  problem.variable_upper_bound(variable_index) - current_iterate.primals[variable_index]);
         }
         // bounds of additional variables (no trust region!)
         for (size_t variable_index: Range(problem.get_number_original_variables(), problem.number_variables)) {
            this->model.lp_.col_lower_[variable_index] = problem.variable_lower_bound(variable_index) - current_iterate.primals[variable_index];
            this->model.lp_.col_upper_[variable_index] = problem.variable_upper_bound(variable_index) - current_iterate.primals[variable_index];
         }
      }

      // constraint bounds
      if (warmstart_information.constraint_bounds_changed || warmstart_information.constraints_changed) {
         for (size_t constraint_index: Range(problem.number_constraints)) {
            this->model.lp_.row_lower_[constraint_index] = problem.constraint_lower_bound(constraint_index) - this->constraints[constraint_index];
            this->model.lp_.row_upper_[constraint_index] = problem.constraint_upper_bound(constraint_index) - this->constraints[constraint_index];
         }
      }

      if (this->print_subproblem) {
         DEBUG << "Linear objective part: "; print_vector(DEBUG, view(this->model.lp_.col_cost_, 0, problem.number_variables));
         DEBUG << "Jacobian:\n";
         DEBUG << "J = "; print_vector(DEBUG, this->model.lp_.a_matrix_.value_);
         DEBUG << "with column start: "; print_vector(DEBUG, this->model.lp_.a_matrix_.start_);
         DEBUG << "and row index: "; print_vector(DEBUG, this->model.lp_.a_matrix_.index_);
         for (size_t variable_index = 0; variable_index < problem.number_variables; variable_index++) {
            DEBUG << "d" << variable_index << " in [" << this->model.lp_.col_lower_[variable_index] << ", " <<
               this->model.lp_.col_upper_[variable_index] << "]\n";
         }
         for (size_t constraint_index = 0; constraint_index < problem.number_constraints; constraint_index++) {
            DEBUG << "linearized c" << constraint_index << " in [" << this->model.lp_.row_lower_[constraint_index] << ", " <<
               this->model.lp_.row_upper_[constraint_index]<< "]\n";
         }
      }
   }

   void HiGHSSolver::solve_subproblem(const OptimizationProblem& problem, Direction& direction) {
      // solve the LP
      HighsStatus return_status = this->highs_solver.passModel(this->model);
      assert(return_status == HighsStatus::kOk);

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
      for (size_t variable_index = 0; variable_index < problem.number_variables; variable_index++) {
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
      for (size_t constraint_index = 0; constraint_index < problem.number_constraints; constraint_index++) {
         direction.multipliers.constraints[constraint_index] = solution.row_dual[constraint_index];
      }
      const HighsInfo& info = this->highs_solver.getInfo();
      direction.subproblem_objective = info.objective_function_value;
   }
} // namespace