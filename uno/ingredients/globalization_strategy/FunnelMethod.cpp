// Copyright (c) 2018-2024 Charlie Vanaret, David Kiessling
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// #include <cmath>
#include "FunnelMethod.hpp"
#include "ProgressMeasures.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

FunnelMethod::FunnelMethod(const Options& options) :
      GlobalizationStrategy(options),
      parameters({
         options.get_double("funnel_ubd"),
         options.get_double("funnel_fact"),
         options.get_double("funnel_delta"),
         options.get_double("funnel_switching_infeasibility_exponent"),
         options.get_double("funnel_kappa_infeasibility_1"),
         options.get_double("funnel_kappa_infeasibility_2"),
         options.get_double("funnel_beta"),
         options.get_double("funnel_gamma"),
         options.get_int("funnel_update_strategy")
      }),
      check_for_current_iterate(options.get_bool("funnel_check_current_point"))
{}

void FunnelMethod::initialize(Statistics& statistics, const Iterate& initial_iterate, const Options& options)
{
   // Add a column of the of the funnel width in the output table
   statistics.add_column("funnel width", Statistics::double_width, options.get_int("statistics_funnel_width_column_order"));

   // set the funnel upper bound
   double upper_bound = std::max(this->parameters.initial_upper_bound,
                                 this->parameters.initial_multiplication * initial_iterate.progress.infeasibility);

   std::cout << "Funnel width: " << upper_bound << std::endl;
   this->funnel_width = upper_bound;
   this->first_iteration_in_solver_phase = true;
   statistics.set("funnel width", this->get_infeasibility_upper_bound());
}

double FunnelMethod::get_infeasibility_upper_bound() const
{
   return this->funnel_width;
}

void FunnelMethod::reset(){}

void FunnelMethod::register_current_progress(const ProgressMeasures& current_progress_measures)
{
   if (this->in_restoration_phase)
   {
    std::cout << "Funnel is reduced after restoration phase" << std::endl;
    this->update_funnel_width_restoration(current_progress_measures.infeasibility);
   }
}

double FunnelMethod::unconstrained_merit_function(const ProgressMeasures& progress)
{
   return progress.objective(1.) + progress.auxiliary;
}

double FunnelMethod::compute_actual_objective_reduction(double current_objective_measure, double /*current_infeasibility_measure*/, double trial_objective_measure)
{
   double actual_reduction = current_objective_measure - trial_objective_measure;
   if (this->protect_actual_reduction_against_roundoff)
   {
      static double machine_epsilon = std::numeric_limits<double>::epsilon();
      actual_reduction += 10. * machine_epsilon * std::abs(current_objective_measure);
   }
   return actual_reduction;
}

bool FunnelMethod::switching_condition(double predicted_reduction, double current_infeasibility, double switching_fraction) const
{
   return predicted_reduction > switching_fraction * std::pow(current_infeasibility, this->parameters.switching_infeasibility_exponent);
}

/*This function checks if the trial iterate infeasibility is inside of the current funnel*/
bool FunnelMethod::is_infeasibility_acceptable_to_funnel(double infeasibility_measure) const
{
   if (infeasibility_measure <= this->parameters.beta*this->funnel_width)
   { // beta is used here
      return true;
   }
   else
   {
      DEBUG << "\t\tNot acceptable to funnel.\n";
      return false;
   }
}

/*This function checks if the funnel sufficient decrease condition is satisfied*/
bool FunnelMethod::is_funnel_sufficient_decrease_satisfied(double infeasibility_measure) const
{
   if (infeasibility_measure <= this->parameters.kappa_infeasibility_1*this->funnel_width)
   { // kappa_1 is used here
      return true;
   }
   else {
      DEBUG << "\t\tFunnel sufficiend decrease not satisfied.\n";
      return false;
   }
}

bool FunnelMethod::is_infeasibility_sufficiently_reduced(const ProgressMeasures& /*current_progress*/, const ProgressMeasures& trial_progress) const
{
   if (this->in_restoration_phase)
   {
      if (this->is_infeasibility_acceptable_to_funnel(trial_progress.infeasibility) && (trial_progress.infeasibility <= this->parameters.kappa_infeasibility_1*this->restoration_entry_infeasibility)){
         return true;
      } else {
         return false;
      }
   }
   else
   {
      return this->is_infeasibility_acceptable_to_funnel(trial_progress.infeasibility);
   }
}

void FunnelMethod::update_funnel_width(double current_infeasibility_measure, double trial_infeasibility_measure)
{
   if (this->parameters.funnel_update_strategy == 1)
   {
      if (trial_infeasibility_measure <= current_infeasibility_measure)
      {
         this->funnel_width = std::max(trial_infeasibility_measure + this->parameters.kappa_infeasibility_2 * (current_infeasibility_measure - trial_infeasibility_measure), this->parameters.kappa_infeasibility_1 *this->funnel_width);
      }
      else
      {
         DEBUG << "Trial infeasibility higher than current infeasibility" << "\n";
         this->funnel_width = this->parameters.kappa_infeasibility_1 *this->funnel_width;
      }
   }
   else if (this->parameters.funnel_update_strategy == 2)
   {
      this->funnel_width = trial_infeasibility_measure + this->parameters.kappa_infeasibility_2 * (this->funnel_width - trial_infeasibility_measure);
      std::cout << "Update 2" << std::endl;
   }
   else if (this->parameters.funnel_update_strategy == 3)
   {
      this->funnel_width = this->parameters.kappa_infeasibility_1 *this->funnel_width;
   }
   std::cout << "Funnel width: " << this->funnel_width << std::endl;
   DEBUG << "\t\tNew funnel parameter is: " << this->funnel_width << "\n";    
}

void FunnelMethod::update_funnel_width_restoration(double current_infeasibility_measure)
{
   this->funnel_width = current_infeasibility_measure + this->parameters.kappa_infeasibility_2 * (this->funnel_width - current_infeasibility_measure);
   std::cout << "Funnel width: " << this->funnel_width << std::endl;
   DEBUG << "\t\tNew funnel parameter is: " << this->funnel_width << "\n"; 
}

bool FunnelMethod::infeasibility_sufficient_reduction(double current_infeasibility, double trial_infeasibility) const {
   return (trial_infeasibility < this->parameters.beta * current_infeasibility);
}

bool FunnelMethod::objective_sufficient_reduction(double current_objective, double trial_objective, double trial_infeasibility) const {
   return (trial_objective <= current_objective - this->parameters.gamma * trial_infeasibility);
}

//! check acceptability wrt current point
bool FunnelMethod::acceptable_wrt_current_iterate(double current_infeasibility, double current_objective, double trial_infeasibility, double trial_objective) {
   return this->objective_sufficient_reduction(current_objective, trial_objective, trial_infeasibility) ||
         this->infeasibility_sufficient_reduction(current_infeasibility, trial_infeasibility);
}

bool FunnelMethod::is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
      const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction, double objective_multiplier)
{
   statistics.set("funnel width", this->get_infeasibility_upper_bound());
   DEBUG << "Current funnel width:\n";
   if (!this->in_restoration_phase && (objective_multiplier == 0.))
   {
      this->restoration_entry_infeasibility = current_progress.infeasibility;
   }
   this->in_restoration_phase = (objective_multiplier == 0.);
   if (this->in_restoration_phase)
   {
      return this->is_feasibility_iterate_acceptable(statistics, current_progress, trial_progress, predicted_reduction);
   }
   else
   {
      return this->is_regular_iterate_acceptable(statistics, current_progress, trial_progress, predicted_reduction);
   }
}

bool FunnelMethod::is_regular_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
         const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction)
{
   bool accept = false;
   std::string scenario;
   bool funnel_acceptable = this->is_infeasibility_acceptable_to_funnel(trial_progress.infeasibility);
   if (funnel_acceptable)
   {
      // in filter and funnel methods, we construct an unconstrained measure by ignoring infeasibility and scaling the objective measure by 1
      const double current_merit = FunnelMethod::unconstrained_merit_function(current_progress);
      const double trial_merit = FunnelMethod::unconstrained_merit_function(trial_progress);
      const double merit_predicted_reduction = FunnelMethod::unconstrained_merit_function(predicted_reduction);
      DEBUG << "Current: (infeasibility, objective + auxiliary) = (" << current_progress.infeasibility << ", " << current_merit << ")\n";
      DEBUG << "Trial:   (infeasibility, objective + auxiliary) = (" << trial_progress.infeasibility << ", " << trial_merit << ")\n";
      DEBUG << "Unconstrained predicted reduction = " << merit_predicted_reduction << '\n';

      // IF check_for_current_iterate == false, then condition always fulfilled, we never check
      // If true, then first part is false, and we always check fur current iterate
      if (!this->check_for_current_iterate || this->acceptable_wrt_current_iterate(current_progress.infeasibility, current_merit, trial_progress.infeasibility, trial_merit))
      {
         // f-type step
         if (this->switching_condition(merit_predicted_reduction, current_progress.infeasibility, this->parameters.delta))
         {
            DEBUG << "\t\tTrial iterate satisfies switching condition ....\n";
            // unconstrained Armijo sufficient decrease condition (predicted reduction should be positive)
            const double objective_actual_reduction = this->compute_actual_objective_reduction(current_merit, current_progress.infeasibility,
                     trial_merit);
            DEBUG << "Unconstrained actual reduction = " << objective_actual_reduction << '\n';
            if (this->armijo_sufficient_decrease(merit_predicted_reduction, objective_actual_reduction))
            {
               DEBUG << "\t\tTrial iterate (f-type) was ACCEPTED by satisfying Armijo condition\n";
               accept = true;
            }
            else
            { // switching condition holds, but not Armijo condition
               DEBUG << "\t\tTrial iterate (f-type) was REJECTED by violating the Armijo condition\n";
            }
            scenario = "f-type Armijo";
         }
         // h-type step
         else if(this->is_funnel_sufficient_decrease_satisfied(trial_progress.infeasibility))
         {
            DEBUG << "\t\tTrial iterate  (h-type) ACCEPTED by violating the switching condition ...\n";
            accept = true; // accept the step and reduce the tr-radius
            DEBUG << "\t\tEntering funnel reduction mechanism\n";
            this->update_funnel_width(current_progress.infeasibility, trial_progress.infeasibility);
            statistics.set("funnel width", this->funnel_width);
            scenario = "h-type";
         }
         else 
         {
            DEBUG << "\t\tTrial iterate REJECTED by violating switching and funnel sufficient decrease condition\n";
            scenario = "none of f/h-type";
         }
      }
      else
      {
         DEBUG << "Trial iterate not acceptable with respect to current point\n";
         scenario = "current point";
      }
   }
   else
   {
      DEBUG << "\t\tTrial iterate REJECTED. Not in funnel\n";
      scenario = "not in funnel";
   }

   statistics.set("status", std::string(accept ? "accepted" : "rejected") + " (" + scenario + ")");
   DEBUG << '\n';
   return accept;
}

bool FunnelMethod::is_feasibility_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
         const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction)
{
   // drop the objective measure and focus on infeasibility and auxiliary terms (barrier, proximal, ...)
   const double current_merit = current_progress.infeasibility + current_progress.auxiliary;
   const double trial_merit = trial_progress.infeasibility + trial_progress.auxiliary;
   const double predicted_merit_reduction = predicted_reduction.infeasibility + predicted_reduction.auxiliary;
   const double actual_merit_reduction = current_merit - trial_merit;
   DEBUG << "Current merit = " << current_merit << '\n';
   DEBUG << "Trial merit = " << trial_merit << '\n';
   DEBUG << "Predicted merit reduction = " << predicted_merit_reduction << '\n';
   DEBUG << "Actual merit reduction = " << actual_merit_reduction << '\n';
   bool accept = false;
   if (this->armijo_sufficient_decrease(predicted_merit_reduction, actual_merit_reduction))
   {
      DEBUG << "Trial iterate (h-type) was accepted by satisfying the Armijo condition\n";
      accept = true;
   }
   else
   {
      DEBUG << "Trial iterate (h-type) was rejected by violating the Armijo condition\n";
   }
   Iterate::number_eval_objective--;
   statistics.set("status", std::string(accept ? "accepted" : "rejected") + " (restoration)");
   return accept;
}
