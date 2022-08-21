// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_AUGMENTEDSYSTEM_H
#define UNO_AUGMENTEDSYSTEM_H

#include <memory>
#include "linear_algebra/SymmetricMatrix.hpp"
#include "optimization/Model.hpp"
#include "solvers/linear/SymmetricIndefiniteLinearSolver.hpp"
#include "tools/Options.hpp"

struct UnstableRegularization : public std::exception {

   [[nodiscard]] const char* what() const noexcept override {
      return "The inertia correction got unstable (delta_w > threshold)";
   }
};

class AugmentedSystem {
public:
   std::unique_ptr<SymmetricMatrix<double>> matrix;
   std::vector<double> rhs{};
   std::vector<double> solution{};

   AugmentedSystem(const std::string& sparse_format, size_t max_dimension, size_t max_number_non_zeros, bool use_regularization,
         const Options& options);
   void factorize_matrix(const Model& model, SymmetricIndefiniteLinearSolver<double>& linear_solver);
   void regularize_matrix(const Model& model, SymmetricIndefiniteLinearSolver<double>& linear_solver, size_t size_primal_block, size_t size_dual_block,
         double dual_regularization_parameter);
   void solve(SymmetricIndefiniteLinearSolver<double>& linear_solver);

protected:
   size_t number_factorizations{0};
   double previous_primal_regularization{0.};
   const double regularization_failure_threshold;
   const double primal_regularization_initial_factor;
   const double dual_regularization_fraction;
   const double primal_regularization_lb;
   const double primal_regularization_decrease_factor;
   const double primal_regularization_fast_increase_factor;
   const double primal_regularization_slow_increase_factor;
};

#endif // UNO_AUGMENTEDSYSTEM_H
