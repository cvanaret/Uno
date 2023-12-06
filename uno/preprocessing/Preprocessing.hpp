// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PREPROCESSING_H
#define UNO_PREPROCESSING_H

#include <vector>
#include "optimization/Model.hpp"
#include "optimization/Iterate.hpp"
#include "solvers/linear/SymmetricIndefiniteLinearSolver.hpp"
#include "tools/Options.hpp"

class Preprocessing {
public:
   static void compute_least_square_multipliers(const Model& variable_index, SymmetricMatrix<double>& matrix, std::vector<double>& rhs,
         SymmetricIndefiniteLinearSolver<double>& linear_solver, Iterate& current_iterate, std::vector<double>& multipliers,
         double multiplier_max_norm);
   static void enforce_linear_constraints(const Options& options, const Model& model, std::vector<double>& x, Multipliers& multipliers);
};

#endif //UNO_PREPROCESSING_H
