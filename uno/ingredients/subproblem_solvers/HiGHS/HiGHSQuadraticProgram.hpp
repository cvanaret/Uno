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
      HiGHSQuadraticProgram(size_t number_variables, size_t number_constraints);

      void initialize_memory(const Subproblem& subproblem) override;
      void build(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) override;
      [[nodiscard]] SolverWorkspace& get_workspace() override;

      HiGHSWorkspace workspace{};
   };
} // namespace

#endif // UNO_HIGHSQUADRATICPROGRAM_H
