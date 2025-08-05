// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BARRIERPARAMETERUPDATESTRATEGY_H
#define UNO_BARRIERPARAMETERUPDATESTRATEGY_H

namespace uno {
   // forward declarations
   class DualResiduals;
   class Iterate;
   class Multipliers;
   class Options;
   class PrimalDualInteriorPointProblem;

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
      [[nodiscard]] bool update_barrier_parameter(const PrimalDualInteriorPointProblem& barrier_problem,
         const Iterate& current_iterate, const Multipliers& current_multipliers, const DualResiduals& residuals);

   protected:
      double barrier_parameter;
      const double tolerance;
      const UpdateParameters parameters;
   };
} // namespace

#endif // UNO_BARRIERPARAMETERUPDATESTRATEGY_H
