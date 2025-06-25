// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BARRIERPARAMETERUPDATESTRATEGY_H
#define UNO_BARRIERPARAMETERUPDATESTRATEGY_H

namespace uno {
   // forward declarations
   class BarrierProblem;
   class DualResiduals;
   class Iterate;
   class Options;

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
      [[nodiscard]] bool update_barrier_parameter(const BarrierProblem& barrier_problem, const Iterate& current_iterate,
         const DualResiduals& residuals);

   protected:
      double barrier_parameter;
      const double dual_tolerance;
      const UpdateParameters parameters;
   };
} // namespace

#endif // UNO_BARRIERPARAMETERUPDATESTRATEGY_H
