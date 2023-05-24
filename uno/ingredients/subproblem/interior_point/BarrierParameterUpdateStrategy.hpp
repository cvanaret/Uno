// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BARRIERPARAMETERUPDATESTRATEGY_H
#define UNO_BARRIERPARAMETERUPDATESTRATEGY_H

#include "reformulation/NonlinearProblem.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Options.hpp"

struct UpdateParameters {
   double k_mu;
   double theta_mu;
   double k_epsilon;
   double update_fraction;
};

class BarrierParameterUpdateStrategy {
public:
   explicit BarrierParameterUpdateStrategy(const Options& options);
   [[nodiscard]] double get_barrier_parameter() const;
   void set_barrier_parameter(double new_barrier_parameter);
   [[nodiscard]] bool update_barrier_parameter(const NonlinearProblem& problem, const Iterate& current_iterate);

protected:
   double barrier_parameter;
   const double tolerance;
   const UpdateParameters parameters;

   [[nodiscard]] static double compute_shifted_complementarity_error(const NonlinearProblem& problem, const Iterate& iterate, double shift_value) ;
};

#endif // UNO_BARRIERPARAMETERUPDATESTRATEGY_H
