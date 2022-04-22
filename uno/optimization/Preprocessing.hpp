#ifndef UNO_PREPROCESSING_H
#define UNO_PREPROCESSING_H

#include <vector>
#include "Model.hpp"
#include "Iterate.hpp"
#include "solvers/linear/LinearSolver.hpp"

class Preprocessing {
public:
   //static void enforce_linear_constraints(const Model& model, Iterate& first_iterate);
   static void compute_least_square_multipliers(const Model& model, SymmetricMatrix& matrix, std::vector<double>& rhs, LinearSolver& solver,
         Iterate& current_iterate, std::vector<double>& multipliers, double multipliers_max_norm = 1e3);
};

#endif //UNO_PREPROCESSING_H
