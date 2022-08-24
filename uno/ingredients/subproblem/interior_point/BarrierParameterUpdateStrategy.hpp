// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BARRIERPARAMETERUPDATESTRATEGY_H
#define UNO_BARRIERPARAMETERUPDATESTRATEGY_H

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
   [[nodiscard]] bool update_barrier_parameter(double primal_dual_error);

protected:
   double barrier_parameter;
   const double tolerance;
   const UpdateParameters parameters;
   const double test{1.};
};

#endif // UNO_BARRIERPARAMETERUPDATESTRATEGY_H
