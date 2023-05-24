// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_RELAXEDPROBLEM_H
#define UNO_RELAXEDPROBLEM_H

#include "NonlinearProblem.hpp"

class RelaxedProblem: public NonlinearProblem {
public:
   RelaxedProblem(const Model& model, size_t number_variables, size_t number_constraints);
   
   [[nodiscard]] virtual double compute_stationarity_error(const Iterate& iterate, Norm residual_norm) const = 0;
   [[nodiscard]] virtual double compute_complementarity_error(const std::vector<double>& primals, const std::vector<double>& constraints,
         const Multipliers& multipliers, Norm residual_norm) const = 0;
};

inline RelaxedProblem::RelaxedProblem(const Model& model, size_t number_variables, size_t number_constraints):
      NonlinearProblem(model, number_variables, number_constraints) {
}

#endif // UNO_RELAXEDPROBLEM_H
