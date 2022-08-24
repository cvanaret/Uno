// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "BarrierParameterUpdateStrategy.hpp"
#include "tools/Logger.hpp"
#include <cmath>

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
   this->barrier_parameter = new_barrier_parameter;
}

bool BarrierParameterUpdateStrategy::update_barrier_parameter(double primal_dual_error) {
   DEBUG << "Max scaled primal-dual error for barrier subproblem is " << primal_dual_error << '\n';

   // update the barrier parameter (Eq. 7 in Ipopt paper)
   bool parameter_updated = false;
   const double tolerance_fraction = this->tolerance / this->parameters.update_fraction;
   while (primal_dual_error <= this->parameters.k_epsilon * this->barrier_parameter && tolerance_fraction < this->barrier_parameter) {
      this->barrier_parameter = std::max(tolerance_fraction, std::min(this->parameters.k_mu * this->barrier_parameter,
            std::pow(this->barrier_parameter, this->parameters.theta_mu)));
      DEBUG << "Barrier parameter mu updated to " << this->barrier_parameter << '\n';
      // the barrier parameter was updated
      parameter_updated = true;
   }
   // the barrier parameter was not updated
   return parameter_updated;
}