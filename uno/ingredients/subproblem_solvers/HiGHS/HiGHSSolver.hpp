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
      explicit HiGHSSolver(const Options& options);

      void initialize_memory(const OptimizationProblem& problem, const HessianModel& hessian_model,
         RegularizationStrategy<double>& regularization_strategy) override;

      void solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
         const Multipliers& current_multipliers, const Vector<double>& initial_point, Direction& direction,
         HessianModel& hessian_model, RegularizationStrategy<double>& regularization_strategy, double trust_region_radius,
         const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] double hessian_quadratic_product(const Vector<double>& vector) const override;

   protected:
      HighsModel model;
      Highs highs_solver;
      std::vector<double> constraints{};
      SparseVector<double> linear_objective{};
      RectangularMatrix<double> constraint_jacobian{};

      const bool print_subproblem;

      void set_up_subproblem(const OptimizationProblem& problem, Iterate& current_iterate, double trust_region_radius,
         const WarmstartInformation& warmstart_information);
      void solve_subproblem(const OptimizationProblem& problem, Direction& direction);
   };
} // namespace

#endif // UNO_HIGHSSOLVER_H