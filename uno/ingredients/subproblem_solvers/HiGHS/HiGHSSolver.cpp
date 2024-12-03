#include <cassert>
#include "HiGHSSolver.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Direction.hpp"
#include "options/Options.hpp"
#include "symbolic/VectorView.hpp"

namespace uno {
   HiGHSSolver::HiGHSSolver(size_t number_variables, size_t number_constraints, size_t number_jacobian_nonzeros, size_t /*number_hessian_nonzeros*/,
         const Options& options): LPSolver(), print_subproblem(options.get_bool("print_subproblem")) {
      this->model.lp_.sense_ = ObjSense::kMinimize;
      this->model.lp_.offset_ = 0.;
      // the linear part of the objective is a dense vector
      this->model.lp_.col_cost_.resize(number_variables);
      // variable bounds
      this->model.lp_.col_lower_.resize(number_variables);
      this->model.lp_.col_upper_.resize(number_variables);
      // constraint bounds
      this->model.lp_.row_lower_.resize(number_constraints);
      this->model.lp_.row_upper_.resize(number_constraints);
      // constraint matrix in CSR format (each row is a constraint gradient)
      this->model.lp_.a_matrix_.format_ = MatrixFormat::kRowwise;
      this->model.lp_.a_matrix_.value_.reserve(number_jacobian_nonzeros);
      this->model.lp_.a_matrix_.index_.reserve(number_jacobian_nonzeros);
      this->model.lp_.a_matrix_.start_.reserve(number_variables + 1);

      this->highs_solver.setOptionValue("output_flag", "false");
   }

   void HiGHSSolver::build_linear_subproblem(size_t number_variables, size_t number_constraints, const std::vector<double>& variables_lower_bounds,
         const std::vector<double>& variables_upper_bounds, const std::vector<double>& constraints_lower_bounds,
         const std::vector<double>& constraints_upper_bounds, const SparseVector<double>& linear_objective,
         const RectangularMatrix<double>& constraint_jacobian) {
      this->model.lp_.num_col_ = static_cast<HighsInt>(number_variables);
      this->model.lp_.num_row_ = static_cast<HighsInt>(number_constraints);

      // variable bounds
      for (size_t variable_index = 0; variable_index < number_variables; variable_index++) {
         this->model.lp_.col_lower_[variable_index] = variables_lower_bounds[variable_index];
         this->model.lp_.col_upper_[variable_index] = variables_upper_bounds[variable_index];
         // reset the linear part of the objective
         this->model.lp_.col_cost_[variable_index] = 0.;
      }

      // linear part of the objective
      for (const auto [variable_index, value]: linear_objective) {
         this->model.lp_.col_cost_[variable_index] = value;
      }

      // constraint bounds
      for (size_t constraint_index = 0; constraint_index < number_constraints; constraint_index++) {
         this->model.lp_.row_lower_[constraint_index] = constraints_lower_bounds[constraint_index];
         this->model.lp_.row_upper_[constraint_index] = constraints_upper_bounds[constraint_index];
      }

      // constraint matrix
      this->model.lp_.a_matrix_.value_.clear();
      this->model.lp_.a_matrix_.index_.clear();
      this->model.lp_.a_matrix_.start_.clear();

      size_t number_nonzeros = 0;
      this->model.lp_.a_matrix_.start_.emplace_back(number_nonzeros);
      for (size_t constraint_index = 0; constraint_index < number_constraints; constraint_index++) {
         for (const auto [variable_index, value]: constraint_jacobian[constraint_index]) {
            this->model.lp_.a_matrix_.value_.emplace_back(value);
            this->model.lp_.a_matrix_.index_.emplace_back(variable_index);
            number_nonzeros++;
         }
         this->model.lp_.a_matrix_.start_.emplace_back(number_nonzeros);
      }

      if (this->print_subproblem) {
         DEBUG << "Linear objective part: "; print_vector(DEBUG, view(this->model.lp_.col_cost_, 0, number_variables));
         DEBUG << "Jacobian:\n";
         DEBUG << "J = "; print_vector(DEBUG, this->model.lp_.a_matrix_.value_);
         DEBUG << "with column start: "; print_vector(DEBUG, this->model.lp_.a_matrix_.start_);
         DEBUG << "and row index: "; print_vector(DEBUG, this->model.lp_.a_matrix_.index_);
         for (size_t variable_index = 0; variable_index < number_variables; variable_index++) {
            DEBUG << "d" << variable_index << " in [" << this->model.lp_.col_lower_[variable_index] << ", " <<
               this->model.lp_.col_upper_[variable_index] << "]\n";
         }
         for (size_t constraint_index = 0; constraint_index < number_constraints; constraint_index++) {
            DEBUG << "linearized c" << constraint_index << " in [" << this->model.lp_.row_lower_[constraint_index] << ", " <<
               this->model.lp_.row_upper_[constraint_index]<< "]\n";
         }
      }
   }

   void HiGHSSolver::solve_subproblem(Direction& direction, size_t number_variables, size_t number_constraints) {
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
      for (size_t variable_index = 0; variable_index < number_variables; variable_index++) {
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
      for (size_t constraint_index = 0; constraint_index < number_constraints; constraint_index++) {
         direction.multipliers.constraints[constraint_index] = solution.row_dual[constraint_index];
      }
      const HighsInfo& info = this->highs_solver.getInfo();
      direction.subproblem_objective = info.objective_function_value;
   }

   void HiGHSSolver::solve_LP(size_t number_variables, size_t number_constraints, const std::vector<double>& variables_lower_bounds,
         const std::vector<double>& variables_upper_bounds, const std::vector<double>& constraints_lower_bounds,
         const std::vector<double>& constraints_upper_bounds, const SparseVector<double>& linear_objective,
         const RectangularMatrix<double>& constraint_jacobian, const Vector<double>& /*initial_point*/, Direction& direction,
         const WarmstartInformation& /*warmstart_information*/) {
      // build the LP in the HiGHS format
      this->build_linear_subproblem(number_variables, number_constraints, variables_lower_bounds, variables_upper_bounds, constraints_lower_bounds,
            constraints_upper_bounds, linear_objective, constraint_jacobian);

      // solve the LP
      this->solve_subproblem(direction, number_variables, number_constraints);
   }
} // namespace