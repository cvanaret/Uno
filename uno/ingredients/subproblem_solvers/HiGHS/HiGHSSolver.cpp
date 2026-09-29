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
      this->highs_solver = Highs_create();
      INFO << "Running HiGHS v" << Highs_version() << '\n';
      Highs_setBoolOptionValue(this->highs_solver, "output_flag", 0);
      // construct an empty HiGHS-native quadratic program so that get_quadratic_program() can be used to
      // build it directly from data (no Subproblem); the full solver instead calls initialize_memory(subproblem)
      this->quadratic_program = std::make_unique<HiGHSQuadraticProgram>();
   }

   HiGHSSolver::~HiGHSSolver() {
      if (this->highs_solver != nullptr) {
         Highs_destroy(this->highs_solver);
      }
   }

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
         DEBUG << "Hessian: " << view(this->quadratic_program->q_value) << '\n';
         DEBUG << "Linear objective part: " << view(this->quadratic_program->col_cost) << '\n';
         DEBUG << "Jacobian: " << view(this->quadratic_program->a_value) << '\n';
         for (size_t variable_index = 0; variable_index < this->quadratic_program->number_variables; variable_index++) {
            DEBUG << "d" << variable_index << " in [" << this->quadratic_program->col_lower[variable_index] << ", " <<
               this->quadratic_program->col_upper[variable_index] << "]\n";
         }
         for (size_t constraint_index = 0; constraint_index < this->quadratic_program->number_constraints; constraint_index++) {
            DEBUG << "linearized c" << constraint_index << " in [" << this->quadratic_program->row_lower[constraint_index] << ", " <<
               this->quadratic_program->row_upper[constraint_index]<< "]\n";
         }
      }
      this->solve_subproblem(direction);
   }

   SolverWorkspace& HiGHSSolver::get_workspace() {
      return *this->quadratic_program;
   }

   // protected member functions

   void HiGHSSolver::solve_subproblem(Direction& direction) {
      HiGHSQuadraticProgram& qp = *this->quadratic_program;
      const HighsInt num_col = static_cast<HighsInt>(qp.number_variables);
      const HighsInt num_row = static_cast<HighsInt>(qp.number_constraints);
      const HighsInt a_num_nz = static_cast<HighsInt>(qp.a_value.size());
      const HighsInt q_num_nz = static_cast<HighsInt>(qp.q_value.size());
      const bool has_hessian = (0 < q_num_nz);

      // pass the model (dense objective + bounds, CSC Jacobian, triangular CSC Hessian).
      // integrality == nullptr: all variables continuous (QP).
      HighsInt return_status = Highs_passModel(this->highs_solver, num_col, num_row, a_num_nz, q_num_nz,
         kHighsMatrixFormatColwise, kHighsHessianFormatTriangular, kHighsObjSenseMinimize, 0.,
         qp.col_cost.data(), qp.col_lower.data(), qp.col_upper.data(),
         qp.row_lower.data(), qp.row_upper.data(),
         qp.a_start.data(), qp.a_index.data(), qp.a_value.data(),
         has_hessian ? qp.q_start.data() : nullptr,
         has_hessian ? qp.q_index.data() : nullptr,
         has_hessian ? qp.q_value.data() : nullptr,
         nullptr);
      if (return_status == kHighsStatusError) {
         throw std::runtime_error("HiGHS could not read the model.");
      }

      DEBUG2 << "Running HiGHS\n";
      return_status = Highs_run(this->highs_solver); // solve
      DEBUG2 << "Ran HiGHS\n";
      DEBUG << "HiGHS status: " << return_status << '\n';

      // if HiGHS could not optimize (e.g. because of indefinite Hessian), return an error
      if (return_status == kHighsStatusError) {
         throw std::runtime_error("HiGHS encountered negative curvature, which it cannot handle. Terminating.");
      }
      const HighsInt model_status = Highs_getModelStatus(this->highs_solver);
      DEBUG << "HiGHS model status: " << model_status << '\n';

      if (model_status == kHighsModelStatusInfeasible) {
         direction.status = SubproblemStatus::INFEASIBLE;
         return;
      }
      else if (model_status == kHighsModelStatusUnbounded) {
         direction.status = SubproblemStatus::UNBOUNDED_PROBLEM;
         return;
      }

      direction.status = SubproblemStatus::OPTIMAL;
      // read the solution into the preallocated workspace buffers
      Highs_getSolution(this->highs_solver, qp.col_value.data(), qp.col_dual.data(), qp.row_value.data(), qp.row_dual.data());
      // read the primal solution and bound dual solution
      for (size_t variable_index = 0; variable_index < qp.number_variables; variable_index++) {
         direction.primals[variable_index] = qp.col_value[variable_index];
         const double bound_multiplier = qp.col_dual[variable_index];
         if (0. < bound_multiplier) {
            direction.multipliers.lower_bounds[variable_index] = bound_multiplier;
         }
         else {
            direction.multipliers.upper_bounds[variable_index] = bound_multiplier;
         }
      }
      // gather the constraint multipliers (the dual-displacement mapping is performed by IQPSolver)
      for (size_t constraint_index = 0; constraint_index < qp.number_constraints; constraint_index++) {
         direction.multipliers.constraints[constraint_index] = qp.row_dual[constraint_index];
      }
      double objective_function_value = 0.;
      Highs_getDoubleInfoValue(this->highs_solver, "objective_function_value", &objective_function_value);
      direction.subproblem_objective = objective_function_value;
   }
} // namespace
