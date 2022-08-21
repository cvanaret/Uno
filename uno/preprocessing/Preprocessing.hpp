// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PREPROCESSING_H
#define UNO_PREPROCESSING_H

#include <vector>
#include "optimization/Model.hpp"
#include "optimization/Iterate.hpp"
#include "solvers/linear/SymmetricIndefiniteLinearSolver.hpp"

class Preprocessing {
public:
   static void compute_least_square_multipliers(const Model& model, SymmetricMatrix<double>& matrix, std::vector<double>& rhs,
         SymmetricIndefiniteLinearSolver<double>& solver, Iterate& current_iterate, std::vector<double>& multipliers,
         double multipliers_max_norm = 1e3);
};

#endif //UNO_PREPROCESSING_H
