// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cassert>
#include "BarrierParameterUpdateStrategy.hpp"
#include "PrimalDualInteriorPointProblem.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Logger.hpp"
#include "options/Options.hpp"

namespace uno {
   BarrierParameterUpdateStrategy::BarrierParameterUpdateStrategy(const Options& options):
      barrier_parameter(options.get_double("barrier_initial_parameter")),
      tolerance(options.get_double("tolerance")),
      parameters({
         options.get_double("barrier_k_mu"),
         options.get_double("barrier_theta_mu"),
         options.get_double("barrier_k_epsilon"),
         options.get_double("barrier_update_fraction")
      }) {
   }

   double BarrierParameterUpdateStrategy::get_barrier_parameter() const {
      return this->barrier_parameter;
   }

   void BarrierParameterUpdateStrategy::set_barrier_parameter(double new_barrier_parameter) {
      assert(0. <= new_barrier_parameter && "The barrier parameter should be positive.");
      this->barrier_parameter = new_barrier_parameter;
   }

   bool BarrierParameterUpdateStrategy::update_barrier_parameter(const PrimalDualInteriorPointProblem& barrier_problem,
         const Iterate& current_iterate, const Multipliers& current_multipliers, const DualResiduals& residuals) {
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
      const double tolerance_fraction = this->tolerance / this->parameters.update_fraction;
      bool parameter_updated = false;
      while (primal_dual_error <= this->parameters.k_epsilon * this->barrier_parameter && tolerance_fraction < this->barrier_parameter) {
         this->barrier_parameter = std::max(tolerance_fraction, std::min(this->parameters.k_mu * this->barrier_parameter,
               std::pow(this->barrier_parameter, this->parameters.theta_mu)));
         DEBUG << "Barrier parameter mu updated to " << this->barrier_parameter << '\n';
         // update complementarity error
         double scaled_complementarity_error = barrier_problem.compute_centrality_error(current_iterate.primals,
            current_multipliers, this->barrier_parameter) / residuals.complementarity_scaling;
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
