// Copyright (c) 2024-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "HiGHSSolver.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
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
      this->evaluation_space.initialize_memory(subproblem);
   }

   void HiGHSSolver::solve(Statistics& statistics, Subproblem& subproblem, double trust_region_radius,
         const Vector<double>& /*initial_point*/, Direction& direction, const WarmstartInformation& warmstart_information) {
      this->set_up_subproblem(statistics, subproblem, trust_region_radius, warmstart_information);
      this->solve_subproblem(subproblem, direction);
   }

   EvaluationSpace& HiGHSSolver::get_evaluation_space() {
      return this->evaluation_space;
   }

   // protected member functions

   void HiGHSSolver::set_up_subproblem(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         const WarmstartInformation& warmstart_information) {
      // evaluate the functions and derivatives
      this->evaluation_space.evaluate_functions(statistics, subproblem, warmstart_information);

      // variable bounds
      if (warmstart_information.variable_bounds_changed) {
         subproblem.set_variables_bounds(this->evaluation_space.model.lp_.col_lower_, this->evaluation_space.model.lp_.col_upper_,
            trust_region_radius);
      }

      // constraint bounds
      if (warmstart_information.constraint_bounds_changed || warmstart_information.new_iterate) {
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