// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HESSIANMODEL_H
#define UNO_HESSIANMODEL_H

#include <memory>
#include <vector>
#include "linear_algebra/SymmetricMatrix.hpp"

namespace uno {
   // forward declarations
   class OptimizationProblem;
   class Statistics;
   template <typename ElementType>
   class Vector;

   class HessianModel {
   public:
      HessianModel(size_t dimension, size_t maximum_number_nonzeros, const std::string& sparse_format, bool use_regularization);
      virtual ~HessianModel();

      SymmetricMatrix<size_t, double> hessian;
      size_t evaluation_count{0};

      virtual void evaluate(Statistics& statistics, const OptimizationProblem& problem, const Vector<double>& primal_variables,
            const Vector<double>& constraint_multipliers) = 0;
   };
} // namespace

#endif // UNO_HESSIANMODEL_H