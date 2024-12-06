#ifndef UNO_HIGHSSOLVER_H
#define UNO_HIGHSSOLVER_H

#include "ingredients/subproblem_solvers/LPSolver.hpp"
#include "Highs.h"
#include "linear_algebra/RectangularMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"

namespace uno {
   // forward declaration
   class Options;

   class HiGHSSolver : public LPSolver {
   public:
      HiGHSSolver(size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros, size_t number_jacobian_nonzeros,
            size_t number_hessian_nonzeros, const Options& options);

      void solve_LP(const OptimizationProblem& problem, Iterate& current_iterate, const Vector<double>& initial_point, Direction& direction,
            double trust_region_radius, const WarmstartInformation& warmstart_information) override;

   protected:
      HighsModel model;
      Highs highs_solver;

      std::vector<double> constraints;
      SparseVector<double> linear_objective;
      RectangularMatrix<double> constraint_jacobian;

      const bool print_subproblem;

      void set_up_subproblem(const OptimizationProblem& problem, Iterate& current_iterate, double trust_region_radius,
            const WarmstartInformation& warmstart_information);
      void solve_subproblem(const OptimizationProblem& problem, Direction& direction);
   };
} // namespace

#endif // UNO_HIGHSSOLVER_H