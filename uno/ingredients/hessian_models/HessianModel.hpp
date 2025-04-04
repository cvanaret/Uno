// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HESSIANMODEL_H
#define UNO_HESSIANMODEL_H

#include <cstddef>
#include <string>

namespace uno {
   // forward declarations
   class Model;
   class Iterate;
   class Options;
   class Statistics;
   template <typename ElementType>
   class Vector;

   class HessianModel {
   public:
      HessianModel() = default;
      virtual ~HessianModel() = default;

      size_t evaluation_count{0};

      [[nodiscard]] virtual bool has_hessian_operator(const Model& model) const = 0;
      [[nodiscard]] virtual bool has_hessian_matrix(const Model& model) const = 0;
      [[nodiscard]] virtual bool has_curvature(const Model& model) const = 0;
      [[nodiscard]] virtual size_t number_nonzeros(const Model& model) const = 0;
      virtual void compute_sparsity(const Model& model, int* row_indices, int* column_indices, int solver_indexing) const = 0;
      [[nodiscard]] virtual bool is_positive_definite() const = 0;

      virtual void initialize(const Model& model) = 0;
      virtual void initialize_statistics(Statistics& statistics, const Options& options) const = 0;
      virtual void evaluate_hessian(Statistics& statistics, const Model& model, const Vector<double>& primal_variables,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* hessian_values) = 0;
      virtual void compute_hessian_vector_product(const Model& model, const double* x, const double* vector,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* result) = 0;
      [[nodiscard]] virtual std::string get_name() const = 0;
   };
} // namespace

#endif // UNO_HESSIANMODEL_H