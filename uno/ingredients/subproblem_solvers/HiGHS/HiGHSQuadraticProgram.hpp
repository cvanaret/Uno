// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HIGHSQUADRATICPROGRAM_H
#define UNO_HIGHSQUADRATICPROGRAM_H

#include "ingredients/subproblem_solvers/QuadraticProgram.hpp"
#include "HiGHSWorkspace.hpp"

namespace uno {
   // forward declarations
   class Evaluations;
   class Statistics;
   class Subproblem;
   class WarmstartInformation;

   // HiGHS-native QuadraticProgram: the data is held directly in a HighsModel (CSC Jacobian + Hessian,
   // dense bounds and objective) inside the workspace. build() fills the model from the Subproblem.
   class HiGHSQuadraticProgram : public QuadraticProgram {
   public:
      HiGHSQuadraticProgram() = default;

      void initialize_memory(const Subproblem& subproblem) override;
      void fill(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) override;
      // data-driven build: dense objective gradient + COO constraint Jacobian + COO Lagrangian Hessian
      // (one triangle; empty for an LP). Converts COO to HiGHS' CSC layout internally.
      void fill(const Vector<double>& linear_objective,
         const Vector<uno_int>& jacobian_row_indices, const Vector<uno_int>& jacobian_column_indices,
         const Vector<double>& jacobian_values,
         const Vector<uno_int>& hessian_row_indices, const Vector<uno_int>& hessian_column_indices,
         const Vector<double>& hessian_values,
         const std::vector<double>& variables_lower_bounds, const std::vector<double>& variables_upper_bounds,
         const std::vector<double>& constraints_lower_bounds, const std::vector<double>& constraints_upper_bounds) override;
      [[nodiscard]] SolverWorkspace& get_workspace() override;

      HiGHSWorkspace workspace{};
   };
} // namespace

#endif // UNO_HIGHSQUADRATICPROGRAM_H
