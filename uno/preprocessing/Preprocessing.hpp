// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PREPROCESSING_H
#define UNO_PREPROCESSING_H

#include <cstddef>
#include <vector>

namespace uno {
   // forward declarations
   template <typename IndexType, typename ElementType>
   class EqualityQPSolver;
   class Iterate;
   class Model;
   class Multipliers;
   class OptimizationProblem;
   class QPSolver;
   class Statistics;
   template <typename IndexType, typename ElementType>
   class SymmetricMatrix;
   template <typename ElementType>
   class Vector;

   class Preprocessing {
   public:
      static void compute_least_square_multipliers(Statistics& statistics, const OptimizationProblem& problem,
         EqualityQPSolver<size_t, double>& equality_QP_solver, Iterate& current_iterate, Multipliers& multipliers,
         double multiplier_max_norm);
      static void enforce_linear_constraints(const Model& model, Vector<double>& primals, Multipliers& multipliers, QPSolver& qp_solver);
   };
} // namespace

#endif //UNO_PREPROCESSING_H
