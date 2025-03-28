// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HESSIANMODEL_H
#define UNO_HESSIANMODEL_H

#include <cstddef>

namespace uno {
   // forward declarations
   class Model;
   class Options;
   class Statistics;
   template <typename IndexType, typename ElementType>
   class SymmetricMatrix;
   template <typename ElementType>
   class Vector;

   class HessianModel {
   public:
      HessianModel() = default;
      virtual ~HessianModel() = default;

      size_t evaluation_count{0};

      [[nodiscard]] virtual size_t compute_number_hessian_nonzeros(const Model& model) const = 0;
      virtual void evaluate(const Model& model, const Vector<double>& primal_variables, double objective_multiplier,
         const Vector<double>& constraint_multipliers, SymmetricMatrix<size_t, double>& hessian) = 0;
   };
} // namespace

#endif // UNO_HESSIANMODEL_H