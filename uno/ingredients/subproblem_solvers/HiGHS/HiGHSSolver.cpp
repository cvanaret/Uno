// Copyright (c) 2024-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <stdexcept>
#include "HiGHSSolver.hpp"
#include "HiGHSQuadraticProgram.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#include "options/Options.hpp"
#include "tools/Infinity.hpp"
#include "tools/Logger.hpp"

namespace uno {
   HiGHSSolver::HiGHSSolver(const Options& options):
         QPSolver(), print_subproblem(options.get_bool("print_subproblem")) {
      INFO << "Running HiGHS v" << highsVersion() << '\n';
      this->highs_solver.setOptionValue("output_flag", "false");
      // construct an empty HiGHS-native quadratic program so that get_quadratic_program() can be used to
      // build it directly from data (no Subproblem); the full solver instead calls initialize_memory(subproblem)
      this->quadratic_program = std::make_unique<HiGHSQuadraticProgram>();
   }

   HiGHSSolver::~HiGHSSolver() = default;

   void HiGHSSolver::initialize_memory(const Subproblem& subproblem) {
      this->quadratic_program->initialize_memory(subproblem);
   }

   QuadraticProgram& HiGHSSolver::get_quadratic_program() {
      return *this->quadratic_program;
   }

   void HiGHSSolver::solve(Statistics& /*statistics*/, const Vector<double>& /*initial_point*/, Direction& direction,
         const WarmstartInformation& /*warmstart_information*/) {
      if (this->print_subproblem) {
         DEBUG << "Subproblem:\n";
         DEBUG << "Hessian: " << view(this->quadratic_program->model.hessian_.value_) << '\n';
         DEBUG << "Linear objective part: " << view(this->quadratic_program->model.lp_.col_cost_) << '\n';
         DEBUG << "Jacobian: " << view(this->quadratic_program->model.lp_.a_matrix_.value_) << '\n';
         for (size_t variable_index = 0; variable_index < this->quadratic_program->number_variables; variable_index++) {
            DEBUG << "d" << variable_index << " in [" << this->quadratic_program->model.lp_.col_lower_[variable_index] << ", " <<
               this->quadratic_program->model.lp_.col_upper_[variable_index] << "]\n";
         }
         for (size_t constraint_index = 0; constraint_index < this->quadratic_program->number_constraints; constraint_index++) {
            DEBUG << "linearized c" << constraint_index << " in [" << this->quadratic_program->model.lp_.row_lower_[constraint_index] << ", " <<
               this->quadratic_program->model.lp_.row_upper_[constraint_index]<< "]\n";
         }
      }
      this->solve_subproblem(direction);
   }

   SolverWorkspace& HiGHSSolver::get_workspace() {
      return *this->quadratic_program;
   }

   // protected member functions

   void HiGHSSolver::solve_subproblem(Direction& direction) {
      // solve the subproblem
      HighsStatus return_status = this->highs_solver.passModel(this->quadratic_program->model);
      if (return_status == HighsStatus::kError) {
         throw std::runtime_error("HiGHS could not read the model.");
      }

      // safeguard: HiGHS' active-set QP solver may stall on convex QPs with a semidefinite Hessian (e.g. the l1 feasibility
      // subproblems, without curvature in the elastic variables). Bound the number of QP iterations generously, so that
      // the solve fails (SubproblemStatus::ERROR) instead of never returning
      const size_t dimension = this->quadratic_program->number_variables + this->quadratic_program->number_constraints;
      this->highs_solver.setOptionValue("qp_iteration_limit", static_cast<HighsInt>(1000 + 10 * dimension));

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

      if (model_status == HighsModelStatus::kInfeasible || model_status == HighsModelStatus::kUnboundedOrInfeasible) {
         direction.status = SubproblemStatus::INFEASIBLE;
         return;
      }
      else if (model_status == HighsModelStatus::kUnbounded) {
         direction.status = SubproblemStatus::UNBOUNDED_PROBLEM;
         return;
      }
      // catch all other errors
      else if (model_status != HighsModelStatus::kOptimal) {
         direction.status = SubproblemStatus::ERROR;
         return;
      }

      direction.status = SubproblemStatus::OPTIMAL;
      const HighsSolution& solution = this->highs_solver.getSolution();
      const HighsBasis& basis = this->highs_solver.getBasis();
      const auto& variables_lower_bounds = this->quadratic_program->model.lp_.col_lower_;
      const auto& variables_upper_bounds = this->quadratic_program->model.lp_.col_upper_;
      // read the primal solution and bound dual solution
      for (size_t variable_index = 0; variable_index < this->quadratic_program->number_variables; variable_index++) {
         // HiGHS' QP solver may report a nonbasic variable at a bound with a value off that bound: place the nonbasic
         // variables on their active bound (consistently with the bound multipliers), and project onto the bounds
         double primal_value = solution.col_value[variable_index];
         if (basis.valid) {
            const HighsBasisStatus variable_status = basis.col_status[variable_index];
            if (variable_status == HighsBasisStatus::kLower && is_finite(variables_lower_bounds[variable_index])) {
               primal_value = variables_lower_bounds[variable_index];
            }
            else if (variable_status == HighsBasisStatus::kUpper && is_finite(variables_upper_bounds[variable_index])) {
               primal_value = variables_upper_bounds[variable_index];
            }
         }
         direction.primals[variable_index] = std::min(std::max(primal_value, variables_lower_bounds[variable_index]),
            variables_upper_bounds[variable_index]);
         const double bound_multiplier = solution.col_dual[variable_index];
         if (0. < bound_multiplier) {
            direction.multipliers.lower_bounds[variable_index] = bound_multiplier;
         }
         else {
            direction.multipliers.upper_bounds[variable_index] = bound_multiplier;
         }
      }
      // gather the constraint multipliers (the dual-displacement mapping is performed by IQPSolver)
      for (size_t constraint_index = 0; constraint_index < this->quadratic_program->number_constraints; constraint_index++) {
         direction.multipliers.constraints[constraint_index] = solution.row_dual[constraint_index];
      }
      const HighsInfo& info = this->highs_solver.getInfo();
      direction.subproblem_objective = info.objective_function_value;
   }
} // namespace
