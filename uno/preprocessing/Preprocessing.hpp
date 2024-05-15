// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PREPROCESSING_H
#define UNO_PREPROCESSING_H

#include <vector>

// forward declarations
class Iterate;
class Model;
class Multipliers;
class QPSolver;
template <typename ElementType>
class SymmetricIndefiniteLinearSolver;
template <typename ElementType>
class SymmetricMatrix;

class Preprocessing {
public:
   static void compute_least_square_multipliers(const Model& model, SymmetricMatrix<double>& matrix, std::vector<double>& rhs,
         SymmetricIndefiniteLinearSolver<double>& linear_solver, Iterate& current_iterate, std::vector<double>& multipliers,
         double multiplier_max_norm);
   [[nodiscard]] static bool enforce_linear_constraints(const Model& model, std::vector<double>& x, Multipliers& multipliers, QPSolver& qp_solver);
};

#endif //UNO_PREPROCESSING_H
