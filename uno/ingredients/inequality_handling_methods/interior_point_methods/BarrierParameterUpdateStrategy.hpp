// Copyright (c) 2018-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BARRIERPARAMETERUPDATESTRATEGY_H
#define UNO_BARRIERPARAMETERUPDATESTRATEGY_H

#include <algorithm>
#include "optimization/DualResiduals.hpp"
#include "optimization/Iterate.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"

namespace uno {
   struct UpdateParameters {
      double k_mu;
      double theta_mu;
      double k_epsilon;
      double update_fraction;
   };

   template <typename BarrierProblem>
   class BarrierParameterUpdateStrategy {
   public:
      explicit BarrierParameterUpdateStrategy(const Options& options);
      [[nodiscard]] double get_barrier_parameter() const;
      void set_barrier_parameter(double new_barrier_parameter);
      [[nodiscard]] bool update_barrier_parameter(const BarrierProblem& barrier_problem,
         const Iterate& current_iterate, const DualResiduals& residuals);

   protected:
      double barrier_parameter;
      const double dual_tolerance;
      const UpdateParameters parameters;
   };

   // class template implementation

   template <typename BarrierProblem>
   BarrierParameterUpdateStrategy<BarrierProblem>::BarrierParameterUpdateStrategy(const Options& options):
      barrier_parameter(options.get_double("barrier_initial_parameter")),
      dual_tolerance(options.get_double("dual_tolerance")),
      parameters({
         options.get_double("barrier_k_mu"),
         options.get_double("barrier_theta_mu"),
         options.get_double("barrier_k_epsilon"),
         options.get_double("barrier_update_fraction")
      }) {
   }

   template <typename BarrierProblem>
   double BarrierParameterUpdateStrategy<BarrierProblem>::get_barrier_parameter() const {
      return this->barrier_parameter;
   }

   template <typename BarrierProblem>
   void BarrierParameterUpdateStrategy<BarrierProblem>::set_barrier_parameter(double new_barrier_parameter) {
      assert(0. <= new_barrier_parameter && "The barrier parameter should be positive.");
      this->barrier_parameter = new_barrier_parameter;
   }

   template <typename BarrierProblem>
   bool BarrierParameterUpdateStrategy<BarrierProblem>::update_barrier_parameter(const BarrierProblem& barrier_problem,
         const Iterate& current_iterate, const DualResiduals& residuals) {
      // primal-dual errors
      const double scaled_stationarity = residuals.stationarity / residuals.stationarity_scaling;
      const double primal_feasibility = (barrier_problem.get_objective_multiplier() == 0.) ? 0. : current_iterate.primal_feasibility;
      double primal_dual_error = std::max({
         scaled_stationarity,
         primal_feasibility,
         residuals.complementarity / residuals.complementarity_scaling
      });
      DEBUG << "Max scaled primal-dual error for barrier subproblem is " << primal_dual_error << '\n';

      // update the barrier parameter (Eq. 7 in IPOPT paper)
      const double tolerance_fraction = this->dual_tolerance / this->parameters.update_fraction;
      bool parameter_updated = false;
      while (primal_dual_error <= this->parameters.k_epsilon * this->barrier_parameter && tolerance_fraction < this->barrier_parameter) {
         this->barrier_parameter = std::max(tolerance_fraction, std::min(this->parameters.k_mu * this->barrier_parameter,
            std::pow(this->barrier_parameter, this->parameters.theta_mu)));
         DEBUG << "Barrier parameter mu updated to " << this->barrier_parameter << '\n';
         // update complementarity error
         double scaled_complementarity_error = barrier_problem.compute_centrality_error(current_iterate.primals,
            current_iterate.multipliers, this->barrier_parameter) / residuals.complementarity_scaling;
         primal_dual_error = std::max({
            scaled_stationarity,
            primal_feasibility,
            scaled_complementarity_error
         });
         DEBUG << "Max scaled primal-dual error for barrier subproblem is " << primal_dual_error << '\n';
         parameter_updated = true;
      }
      return parameter_updated;
   }
} // namespace

#endif // UNO_BARRIERPARAMETERUPDATESTRATEGY_H
