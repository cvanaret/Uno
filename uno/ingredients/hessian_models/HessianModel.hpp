// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HESSIANMODEL_H
#define UNO_HESSIANMODEL_H

#include <memory>
#include <vector>
#include "linear_algebra/SymmetricMatrix.hpp"
#include "solvers/DirectSymmetricIndefiniteLinearSolver.hpp"

namespace uno {
   // forward declarations
   /*
   template <typename IndexType, typename NumericalType>
   class DirectSymmetricIndefiniteLinearSolver;
    */
   class OptimizationProblem;
   class Options;
   class Statistics;
   template <typename ElementType>
   class Vector;

   class HessianModel {
   public:
      HessianModel(size_t dimension, size_t maximum_number_nonzeros, const std::string& sparse_format, bool use_regularization);
      virtual ~HessianModel();

      SymmetricMatrix<size_t, double> hessian;
      size_t evaluation_count{0};

      virtual void evaluate(const OptimizationProblem& problem, const Vector<double>& primals, const Vector<double>& constraint_multipliers,
            SymmetricMatrix<size_t, double>& hessian) = 0;
   };

   // Exact Hessian
   class ExactHessian : public HessianModel {
   public:
      ExactHessian(size_t dimension, size_t maximum_number_nonzeros, const Options& options);

      void evaluate(const OptimizationProblem& problem, const Vector<double>& primals, const Vector<double>& constraint_multipliers,
            SymmetricMatrix<size_t, double>& hessian) override;
   };

   // Hessian with convexification (inertia correction)
   class ConvexifiedHessian : public HessianModel {
   public:
      ConvexifiedHessian(size_t dimension, size_t maximum_number_nonzeros, const Options& options);

      void evaluate(const OptimizationProblem& problem, const Vector<double>& primals, const Vector<double>& constraint_multipliers,
            SymmetricMatrix<size_t, double>& hessian) override;

   protected:
      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<size_t, double>> linear_solver; /*!< Solver that computes the inertia */
      const double regularization_initial_value{};
      const double regularization_increase_factor{};

      void regularize(SymmetricMatrix<size_t, double>& hessian, size_t number_original_variables);
   };

   // zero Hessian
   class ZeroHessian : public HessianModel {
   public:
      ZeroHessian(size_t dimension, size_t maximum_number_nonzeros, const Options& options);

      void evaluate(const OptimizationProblem& problem, const Vector<double>& primals, const Vector<double>& constraint_multipliers,
            SymmetricMatrix<size_t, double>& hessian) override;
   };
} // namespace

#endif // UNO_HESSIANMODEL_H