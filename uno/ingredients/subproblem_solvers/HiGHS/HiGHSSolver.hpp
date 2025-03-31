#ifndef UNO_HIGHSSOLVER_H
#define UNO_HIGHSSOLVER_H

#include "ingredients/subproblem_solvers/LPSolver.hpp"
#include "Highs.h"
#include "linear_algebra/RectangularMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declaration
   class Options;

   class HiGHSSolver : public LPSolver {
   public:
      HiGHSSolver(size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros, size_t number_jacobian_nonzeros,
            size_t number_hessian_nonzeros, const Options& options);

      [[nodiscard]] SubproblemStatus solve_LP(Statistics& statistics, LagrangeNewtonSubproblem& subproblem,
         const Vector<double>& initial_point, Vector<double>& direction_primals, Multipliers& direction_multipliers, double& subproblem_objective,
         const WarmstartInformation& warmstart_information) override;

   protected:
      HighsModel model;
      Highs highs_solver;

      Vector<double> constraints;
      SparseVector<double> linear_objective;
      RectangularMatrix<double> constraint_jacobian;

      const bool print_subproblem;

      void set_up_subproblem(LagrangeNewtonSubproblem& subproblem, const WarmstartInformation& warmstart_information);
      [[nodiscard]] SubproblemStatus solve_subproblem(LagrangeNewtonSubproblem& subproblem, Vector<double>& direction_primals,
         Multipliers& direction_multipliers, double& subproblem_objective);
   };
} // namespace

#endif // UNO_HIGHSSOLVER_H