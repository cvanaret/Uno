// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <memory>
#include "HessianModel.hpp"

namespace uno {
   // forward declarations
   template <typename IndexType, typename NumericalType>
   class DirectSymmetricIndefiniteLinearSolver;
   class Options;

   // Hessian with convexification (inertia correction)
   class ConvexifiedHessian : public HessianModel {
   public:
      ConvexifiedHessian(size_t dimension, size_t maximum_number_nonzeros, const Options& options);

      void evaluate(Statistics& statistics, const OptimizationProblem& problem, const Vector<double>& primal_variables,
            const Vector<double>& constraint_multipliers) override;

   protected:
      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<size_t, double>> linear_solver; /*!< Solver that computes the inertia */
      const double regularization_initial_value{};
      const double regularization_increase_factor{};
      const double regularization_failure_threshold{};

      void regularize(Statistics& statistics, SymmetricMatrix<size_t, double>& hessian, size_t number_original_variables);
   };
} // namespace