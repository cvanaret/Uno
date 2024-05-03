// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_RELAXEDPROBLEM_H
#define UNO_RELAXEDPROBLEM_H

#include "OptimizationProblem.hpp"

class RelaxedProblem: public OptimizationProblem {
public:
   RelaxedProblem(const Model& model, size_t number_variables, size_t number_constraints);
   ~RelaxedProblem() override = default;
   
   [[nodiscard]] virtual double stationarity_error(const LagrangianGradient<double>& lagrangian_gradient, Norm residual_norm) const = 0;
   [[nodiscard]] virtual double complementarity_error(const std::vector<double>& primals, const std::vector<double>& constraints,
         const Multipliers& multipliers, Norm residual_norm) const = 0;
};

inline RelaxedProblem::RelaxedProblem(const Model& model, size_t number_variables, size_t number_constraints):
      OptimizationProblem(model, number_variables, number_constraints) {
}

#endif // UNO_RELAXEDPROBLEM_H
