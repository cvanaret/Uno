// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HIGHSSOLVER_H
#define UNO_HIGHSSOLVER_H

#include <memory>
#include <highs/Highs.h>
#include "../QPSolver.hpp"

namespace uno {
   // forward declarations
   class HiGHSQuadraticProgram;
   class Options;
   class QuadraticProgram;

   class HiGHSSolver : public QPSolver {
   public:
      explicit HiGHSSolver(const Options& options);
      ~HiGHSSolver() override;

      void initialize_memory(const Subproblem& subproblem) override;

      [[nodiscard]] QuadraticProgram& get_quadratic_program() override;

      void solve(Statistics& statistics, const Vector<double>& initial_point, Direction& direction,
         const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] SolverWorkspace& get_workspace() override;

   protected:
      // HiGHS-native quadratic program (built by IQPSolver before each solve)
      std::unique_ptr<HiGHSQuadraticProgram> quadratic_program{};
      Highs highs_solver;

      const bool print_subproblem;

      void solve_subproblem(Direction& direction);
   };
} // namespace

#endif // UNO_HIGHSSOLVER_H
