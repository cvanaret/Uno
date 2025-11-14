// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HIGHSSOLVER_H
#define UNO_HIGHSSOLVER_H

#include "ingredients/subproblem_solvers/QPSolver.hpp"
#include "Highs.h"
#include "HiGHSEvaluationSpace.hpp"

namespace uno {
   // forward declaration
   class Options;

   class HiGHSSolver : public QPSolver {
   public:
      explicit HiGHSSolver(const Options& options);

      void initialize_memory(const Subproblem& subproblem) override;

      void solve(Statistics& statistics, Subproblem& subproblem, double trust_region_radius, const Vector<double>& initial_point,
         Direction& direction, const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] EvaluationSpace& get_evaluation_space() override;

   protected:
      Highs highs_solver;
      HiGHSEvaluationSpace evaluation_space;

      const bool print_subproblem;

      void set_up_subproblem(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         const WarmstartInformation& warmstart_information);
      void solve_subproblem(const Subproblem& subproblem, Direction& direction);
   };
} // namespace

#endif // UNO_HIGHSSOLVER_H