// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HESSIANMODEL_H
#define UNO_HESSIANMODEL_H

#include <cstddef>
#include <string>
#include <string_view>

namespace uno {
   // forward declarations
   class Statistics;
   template <typename ElementType>
   class Vector;

   class HessianModel {
   public:
      explicit HessianModel(const std::string_view name): name(name) {
      }
      virtual ~HessianModel() = default;

      size_t evaluation_count{0};
      const std::string name;

      [[nodiscard]] virtual bool has_hessian_operator() const = 0;
      [[nodiscard]] virtual bool has_hessian_matrix() const = 0;
      [[nodiscard]] virtual bool has_curvature() const = 0;
      [[nodiscard]] virtual size_t number_nonzeros() const = 0;
      virtual void compute_sparsity(int* row_indices, int* column_indices, int solver_indexing) const = 0;
      [[nodiscard]] virtual bool is_positive_definite() const = 0;

      virtual void evaluate_hessian(Statistics& statistics, const Vector<double>& primal_variables,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* hessian_values) = 0;
      virtual void compute_hessian_vector_product(const double* x, const double* vector,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* result) = 0;
   };
} // namespace

#endif // UNO_HESSIANMODEL_H