// Copyright (c) 2024-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "HiGHSSolver.hpp"
#include "HiGHSQuadraticProgram.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"

namespace uno {
   HiGHSSolver::HiGHSSolver(const Options& options):
         QPSolver(), print_subproblem(options.get_bool("print_subproblem")) {
      this->highs_solver.setOptionValue("output_flag", "false");
   }

   HiGHSSolver::~HiGHSSolver() = default;

   void HiGHSSolver::initialize_memory(const Subproblem& subproblem) {
      this->quadratic_program = std::make_unique<HiGHSQuadraticProgram>(subproblem.number_variables,
         subproblem.number_constraints);
      this->quadratic_program->initialize_memory(subproblem);
   }

   QuadraticProgram& HiGHSSolver::get_quadratic_program() {
      return *this->quadratic_program;
   }

   void HiGHSSolver::solve(Statistics& /*statistics*/, const Vector<double>& /*initial_point*/, Direction& direction,
         const WarmstartInformation& /*warmstart_information*/) {
      if (this->print_subproblem) {
         const HiGHSQuadraticProgram& quadratic_program = *this->quadratic_program;
         DEBUG << "Subproblem:\n";
         DEBUG << "Hessian: "; print_vector(DEBUG, quadratic_program.workspace.model.hessian_.value_);
         DEBUG << "Linear objective part: "; print_vector(DEBUG, quadratic_program.workspace.model.lp_.col_cost_);
         DEBUG << "Jacobian: "; print_vector(DEBUG, quadratic_program.workspace.model.lp_.a_matrix_.value_);
         for (size_t variable_index = 0; variable_index < quadratic_program.number_variables; variable_index++) {
            DEBUG << "d" << variable_index << " in [" << quadratic_program.workspace.model.lp_.col_lower_[variable_index] << ", " <<
               quadratic_program.workspace.model.lp_.col_upper_[variable_index] << "]\n";
         }
         for (size_t constraint_index = 0; constraint_index < quadratic_program.number_constraints; constraint_index++) {
            DEBUG << "linearized c" << constraint_index << " in [" << quadratic_program.workspace.model.lp_.row_lower_[constraint_index] << ", " <<
               quadratic_program.workspace.model.lp_.row_upper_[constraint_index]<< "]\n";
         }
      }
      this->solve_subproblem(direction);
   }

   SolverWorkspace& HiGHSSolver::get_workspace() {
      return this->quadratic_program->get_workspace();
   }

   // protected member functions

   void HiGHSSolver::solve_subproblem(Direction& direction) {
      HiGHSQuadraticProgram& quadratic_program = *this->quadratic_program;

      // solve the subproblem
      HighsStatus return_status = this->highs_solver.passModel(quadratic_program.workspace.model);
      if (return_status == HighsStatus::kError) {
         throw std::runtime_error("HiGHS could not read the model.");
      }

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
      for (size_t variable_index = 0; variable_index < quadratic_program.number_variables; variable_index++) {
         direction.primals[variable_index] = solution.col_value[variable_index];
         const double bound_multiplier = solution.col_dual[variable_index];
         if (0. < bound_multiplier) {
            direction.multipliers.lower_bounds[variable_index] = bound_multiplier;
         }
         else {
            direction.multipliers.upper_bounds[variable_index] = bound_multiplier;
         }
      }
      // gather the constraint multipliers (the dual-displacement mapping is performed by IQPSolver)
      for (size_t constraint_index = 0; constraint_index < quadratic_program.number_constraints; constraint_index++) {
         direction.multipliers.constraints[constraint_index] = solution.row_dual[constraint_index];
      }
      const HighsInfo& info = this->highs_solver.getInfo();
      direction.subproblem_objective = info.objective_function_value;
   }
} // namespace
