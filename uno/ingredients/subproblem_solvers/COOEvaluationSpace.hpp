// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_COOEVALUATIONSPACE_H
#define UNO_COOEVALUATIONSPACE_H

#include <cstddef>
#include <vector>
#include "linear_algebra/Vector.hpp"
#include "optimization/EvaluationSpace.hpp"

namespace uno {
   class COOEvaluationSpace: public EvaluationSpace {
   public:
      COOEvaluationSpace() = default;
      ~COOEvaluationSpace() override = default;

      void initialize_hessian(const Subproblem& subproblem);
      void initialize_augmented_system(const Subproblem& subproblem);

      void evaluate_constraint_jacobian(const OptimizationProblem& problem, Iterate& iterate) override;
      void compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const override;
      void compute_constraint_jacobian_transposed_vector_product(const Vector<double>& vector,
         Vector<double>& result) const override;
      [[nodiscard]] double compute_hessian_quadratic_product(const Vector<double>& vector) const override;

      void set_up_linear_system(Statistics& statistics, const Subproblem& subproblem, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         const WarmstartInformation& warmstart_information);

      Vector<double> objective_gradient{}; /*!< Sparse Jacobian of the objective */
      Vector<double> constraints{}; /*!< Constraint values (size \f$m)\f$ */

      // Jacobian
      size_t number_jacobian_nonzeros{};
      std::vector<int> jacobian_row_indices{};
      std::vector<int> jacobian_column_indices{};

      // symmetric matrix (Hessian or augmented system)
      size_t number_hessian_nonzeros{};
      size_t number_matrix_nonzeros{};
      std::vector<int> matrix_row_indices{};
      std::vector<int> matrix_column_indices{};
      Vector<double> matrix_values;
      Vector<double> rhs{};
      Vector<double> solution{};
      bool analysis_performed{false};
   };
} // namespace

#endif // UNO_COOEVALUATIONSPACE_H