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

      void solve_LP(size_t number_variables, size_t number_constraints, const std::vector<double>& variables_lower_bounds,
            const std::vector<double>& variables_upper_bounds, const std::vector<double>& constraints_lower_bounds,
            const std::vector<double>& constraints_upper_bounds, const SparseVector<double>& linear_objective,
            const RectangularMatrix<double>& constraint_jacobian, const Vector<double>& initial_point, Direction& direction,
            const WarmstartInformation& warmstart_information) override;

   protected:
      HighsModel model;
      Highs highs_solver;
      const bool print_subproblem;

      void build_linear_subproblem(size_t number_variables, size_t number_constraints, const std::vector<double>& variables_lower_bounds,
            const std::vector<double>& variables_upper_bounds, const std::vector<double>& constraints_lower_bounds,
            const std::vector<double>& constraints_upper_bounds, const SparseVector<double>& linear_objective,
            const RectangularMatrix<double>& constraint_jacobian);
      void solve_subproblem(Direction& direction, size_t number_variables, size_t number_constraints);
   };
} // namespace

#endif // UNO_HIGHSSOLVER_H