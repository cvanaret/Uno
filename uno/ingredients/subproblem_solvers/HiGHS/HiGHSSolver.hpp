// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HIGHSSOLVER_H
#define UNO_HIGHSSOLVER_H

#include "ingredients/subproblem_solvers/LPSolver.hpp"
#include "Highs.h"
#include "linear_algebra/COOFormat.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "linear_algebra/SparseSymmetricMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"

namespace uno {
   // forward declaration
   class Options;

   class HiGHSSolver : public LPSolver {
   public:
      explicit HiGHSSolver(const Options& options);

      void initialize_memory(const OptimizationProblem& problem, const HessianModel& hessian_model,
         const RegularizationStrategy<double>& regularization_strategy) override;

      void solve(Statistics& statistics, Subproblem& subproblem, const Vector<double>& initial_point,
         Direction& direction, const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] double hessian_quadratic_product(const Vector<double>& vector) const override;

   protected:
      HighsModel model;
      Highs highs_solver;
      std::vector<double> constraints{};
      SparseVector<double> linear_objective{};
      RectangularMatrix<double> constraint_jacobian{};
      SparseSymmetricMatrix<COOFormat<size_t, double>> hessian{};

      const bool print_subproblem;

      void set_up_subproblem(Statistics& statistics, const Subproblem& subproblem, const WarmstartInformation& warmstart_information);
      void solve_subproblem(const Subproblem& subproblem, Direction& direction);
   };
} // namespace

#endif // UNO_HIGHSSOLVER_H