#include <cmath>
#include <assert.h>
#include "TrustRegion.hpp"
#include "linear_algebra/Vector.hpp"
#include "tools/Logger.hpp"

TrustRegion::TrustRegion(ConstraintRelaxationStrategy& constraint_relaxation_strategy, double initial_radius, int max_iterations) :
      GlobalizationMechanism(constraint_relaxation_strategy, max_iterations), radius(initial_radius) {
}

Iterate TrustRegion::initialize(Statistics& statistics, const Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
   statistics.add_column("TR radius", Statistics::double_width, 30);
   // generate the initial point
   Iterate first_iterate = this->relaxation_strategy.initialize(statistics, problem, x, multipliers);

   // preallocate trial_iterate
   this->trial_primals_.resize(first_iterate.x.size());
   this->trial_duals_.resize(multipliers.constraints.size());

   return first_iterate;
}

std::tuple<Iterate, double, double> TrustRegion::compute_acceptable_iterate(Statistics& statistics, const Problem& problem, Iterate&
current_iterate) {
   this->number_iterations = 0;

   while (!this->termination()) {
      try {
         assert (0 < this->radius);
         this->number_iterations++;
         this->print_iteration();

         /* generate the subproblem */
         this->relaxation_strategy.generate_subproblem(problem, current_iterate, this->radius);

         /* compute the direction within the trust region */
         Direction direction = this->relaxation_strategy.compute_feasible_direction(statistics, problem, current_iterate);
         /* set bound multipliers of active trust region to 0 */
         TrustRegion::rectify_active_set(direction, this->radius);

         // assemble the trial iterate
         Iterate trial_iterate = this->assemble_trial_iterate(current_iterate, direction, 1.);

         // check whether the trial step is accepted
         if (this->relaxation_strategy.is_acceptable(statistics, problem, current_iterate, trial_iterate, direction, 1.)) {
            this->add_statistics(statistics, direction);

            // increase the radius if trust region is active
            if (direction.norm >= this->radius - this->activity_tolerance) {
               this->radius *= this->increase_factor;
            }

            // let the subproblem know the accepted iterate
            this->relaxation_strategy.register_accepted_iterate(trial_iterate);
            return std::make_tuple(std::move(trial_iterate), direction.norm, direction.objective_multiplier);
         }
         else {
            /* if the step is rejected, decrease the radius */
            this->radius = std::min(this->radius, direction.norm) / this->decrease_factor;
         }
      }
      catch (const NumericalError& e) {
         GlobalizationMechanism::print_warning(e.what());
         /* if an evaluation error occurs, decrease the radius */
         this->radius /= this->decrease_factor;
      }
   }
   if (this->max_iterations < this->number_iterations) {
      throw std::runtime_error("Trust-region iteration limit reached");
   } /* radius gets too small */
   else if (this->radius < this->min_radius) { /* 1e-16: something like machine precision */
      throw std::runtime_error("Trust-region radius became too small");
   }
   else {
      throw std::runtime_error("Trust-region failed with an unexpected error");
   }
}

void TrustRegion::add_statistics(Statistics& statistics, const Direction& direction) {
   statistics.add_statistic("minor", this->number_iterations);
   statistics.add_statistic("TR radius", this->radius);
   statistics.add_statistic("step norm", direction.norm);
}

void TrustRegion::rectify_active_set(Direction& direction, double radius) {
   assert (0 < radius);
   /* update active set and set multipliers for bound constraints active at trust region to 0 */
   for (auto it = direction.active_set.bounds.at_lower_bound.begin(); it != direction.active_set.bounds.at_lower_bound.end();) {
      int i = *it;
      if (direction.x[i] == -radius) {
         it = direction.active_set.bounds.at_lower_bound.erase(it);
         direction.multipliers.lower_bounds[i] = 0.;
      }
      else {
         ++it;
      }
   }
   for (auto it = direction.active_set.bounds.at_upper_bound.begin(); it != direction.active_set.bounds.at_upper_bound.end();) {
      int i = *it;
      if (direction.x[i] == radius) {
         it = direction.active_set.bounds.at_upper_bound.erase(it);
         direction.multipliers.upper_bounds[i] = 0.;
      }
      else {
         ++it;
      }
   }
}

bool TrustRegion::termination() {
   return (this->max_iterations < this->number_iterations || this->radius < this->min_radius);
}

void TrustRegion::print_iteration() {
   DEBUG << "\n\tTRUST REGION iteration " << this->number_iterations << ", radius " << this->radius << "\n";
}
