#ifndef UNO_HIGHSSOLVER_H
#define UNO_HIGHSSOLVER_H

#include "solvers/LPSolver.hpp"
#include "Highs.h"

namespace uno {
   // forward declaration
   class Options;

   class HiGHSSolver : public LPSolver {
   public:
      HiGHSSolver(size_t number_variables, size_t number_constraints, size_t number_jacobian_nonzeros, size_t number_hessian_nonzeros,
            const Options& options);

      void solve_LP(const LagrangeNewtonSubproblem& subproblem, const Vector<double>& initial_point, Direction& direction,
            const WarmstartInformation& warmstart_information) override;

   protected:
      HighsModel model;
      Highs highs_solver;
      const bool print_subproblem;

      void build_linear_subproblem(const LagrangeNewtonSubproblem& subproblem, const WarmstartInformation& warmstart_information);
      void solve_subproblem(const LagrangeNewtonSubproblem& subproblem, Direction& direction);
   };
} // namespace

#endif // UNO_HIGHSSOLVER_H