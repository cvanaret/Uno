// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_AUGMENTEDSYSTEM_H
#define UNO_AUGMENTEDSYSTEM_H

#include <memory>
#include "linear_algebra/SymmetricMatrix.hpp"
#include "optimization/Model.hpp"
#include "solvers/linear/LinearSolver.hpp"
#include "tools/Options.hpp"

struct UnstableRegularization : public std::exception {

   [[nodiscard]] const char* what() const noexcept override {
      return "The inertia correction got unstable (delta_w > threshold)";
   }
};

class AugmentedSystem {
public:
   std::unique_ptr<SymmetricMatrix> matrix;
   std::vector<double> rhs{};
   std::vector<double> solution{};

   AugmentedSystem(const std::string& sparse_format, size_t max_dimension, size_t max_number_non_zeros, bool use_regularization,
         const Options& options);
   void factorize_matrix(const Model& model, LinearSolver& linear_solver);
   void regularize_matrix(const Model& model, LinearSolver& linear_solver, size_t size_first_block, size_t size_second_block,
         double constraint_regularization_parameter);
   void solve(LinearSolver& linear_solver);

protected:
   size_t number_factorizations{0};
   double previous_regularization_first_block{0.};
   const double regularization_failure_threshold;
   const double regularization_first_block_initial_factor;
   const double regularization_second_block_fraction;
   const double regularization_first_block_lb;
   const double regularization_first_block_decrease_factor;
   const double regularization_first_block_fast_increase_factor;
   const double regularization_first_block_slow_increase_factor;
};

#endif // UNO_AUGMENTEDSYSTEM_H
