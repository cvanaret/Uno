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

      void initialize_statistics(Statistics& statistics, const Options& options) const override;
      void evaluate_hessian(const Model& model, const Vector<double>& primal_variables, double objective_multiplier,
         const Vector<double>& constraint_multipliers, SymmetricMatrix<size_t, double>& hessian) override;
      void compute_hessian_vector_product(const Model& model, const Vector<double>& vector, double objective_multiplier,
         const Vector<double>& constraint_multipliers, Vector<double>& result) override;

   protected:
      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<size_t, double>> linear_solver; /*!< Solver that computes the inertia */
      double regularization_factor{0.};
      const double regularization_initial_value{};
      const double regularization_increase_factor{};
      const double regularization_failure_threshold{};

      void regularize(Statistics& statistics, SymmetricMatrix<size_t, double>& hessian, size_t number_original_variables);
   };
} // namespace