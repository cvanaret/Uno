// Copyright (c) 2018-2024 Charlie Vanaret, David Kiessling
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// #include <cmath>
#include "FunnelMethod.hpp"
#include "ProgressMeasures.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Logger.hpp"
#include "tools/Options.hpp"
#include "tools/Statistics.hpp"

FunnelMethod::FunnelMethod(const Options& options) :
      GlobalizationStrategy(options),
      parameters({
         options.get_double("funnel_kappa_initial_upper_bound"),
         options.get_double("funnel_kappa_initial_multiplication"),
         options.get_double("funnel_delta"),
         options.get_double("funnel_ubd"),
         options.get_double("funnel_fact"),
         options.get_double("funnel_switching_infeasibility_exponent"),
         options.get_double("funnel_kappa_infeasibility_1"),
         options.get_double("funnel_kappa_infeasibility_2"),
         options.get_double("funnel_beta")
      })
      {
}

void FunnelMethod::initialize(Statistics& statistics, const Iterate& initial_iterate, const Options& options) {
    // Add a column of the of the funnel width in the output table
    statistics.add_column("funnel width", Statistics::double_width, options.get_int("statistics_funnel_width_column_order"));

   // set the funnel upper bound
   double upper_bound = std::max(this->parameters.kappa_initial_upper_bound,
                                 this->parameters.kappa_initial_multiplication * initial_iterate.progress.infeasibility);

   this->funnel_width = upper_bound;
   this->first_iteration_in_solver_phase = true;
   statistics.set("funnel width", this->get_funnel_width());
}

void FunnelMethod::reset(){}

void FunnelMethod::register_current_progress(const ProgressMeasures& current_progress_measures) {
   if (!this->in_restoration_phase){
    std::cout << "Funnel is reduced after restoration phase" << std::endl;
    this->update_funnel_width_restoration(current_progress_measures.infeasibility);
   }
}

double FunnelMethod::get_infeasibility_upper_bound() const {
    return this->funnel_width;
}

double FunnelMethod::compute_actual_objective_reduction(double current_objective_measure, double /*current_infeasibility_measure*/, double trial_objective_measure) {
   double actual_reduction = current_objective_measure - trial_objective_measure;
   if (this->protect_actual_reduction_against_roundoff) {
      static double machine_epsilon = std::numeric_limits<double>::epsilon();
      actual_reduction += 10. * machine_epsilon * std::abs(current_objective_measure);
   }
   return actual_reduction;
}

bool FunnelMethod::switching_condition(double predicted_reduction, double current_infeasibility, double switching_fraction) const {
   return predicted_reduction > switching_fraction * std::pow(current_infeasibility, this->parameters.switching_infeasibility_exponent);
}

/*This function checks if the trial iterate infeasibility is inside of the current funnel*/
bool FunnelMethod::is_infeasibility_acceptable_to_funnel(double infeasibility_measure) const {
   if (infeasibility_measure <= this->parameters.beta*this->funnel_width){ // beta is used here
      return true;
   }
   else {
      DEBUG << "\t\tNot acceptable to funnel.\n";
      return false;
   }
}

/*This function checks if the funnel sufficient decrease condition is satisfied*/
bool FunnelMethod::is_funnel_sufficient_decrease_satisfied(double infeasibility_measure) const {
   if (infeasibility_measure <= this->parameters.kappa_infeasibility_1*this->funnel_width){ // kappa_1 is used here
      return true;
   }
   else {
      DEBUG << "\t\tFunnel sufficiend decrease not satisfied.\n";
      return false;
   }
}

bool FunnelMethod::is_infeasibility_sufficiently_reduced(const ProgressMeasures& /*current_progress*/, const ProgressMeasures& trial_progress) const {
   if (this->in_restoration_phase){
        if (this->is_infeasibility_acceptable_to_funnel(trial_progress.infeasibility) && (trial_progress.infeasibility <= this->parameters.kappa_infeasibility_1*this->restoration_entry_infeasibility)){
            return true;
        } else {
            return false;
        }
   } else {
        return this->is_infeasibility_acceptable_to_funnel(trial_progress.infeasibility);
   }
}

void FunnelMethod::update_funnel_width(double current_infeasibility_measure, double trial_infeasibility_measure) {

    this->funnel_width = std::max(this->parameters.kappa_infeasibility_1 *this->funnel_width, 
      trial_infeasibility_measure + this->parameters.kappa_infeasibility_2 * (current_infeasibility_measure - trial_infeasibility_measure));
   DEBUG << "\t\tNew funnel parameter is: " << this->funnel_width << "\n"; 
   
}

void FunnelMethod::update_funnel_width_restoration(double current_infeasibility_measure) {
   this->funnel_width = std::min(current_infeasibility_measure + this->parameters.kappa_infeasibility_2 * (this->funnel_width - current_infeasibility_measure), this->parameters.beta*this->funnel_width);
   DEBUG << "\t\tNew funnel parameter is: " << this->funnel_width << "\n"; 
}

double FunnelMethod::get_funnel_width(){
   return this->funnel_width;
}

double FunnelMethod::unconstrained_merit_function(const ProgressMeasures& progress) {
   return progress.objective(1.) + progress.auxiliary;
}

bool FunnelMethod::is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
         const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction, double objective_multiplier) {
   this->in_restoration_phase = (objective_multiplier == 0.);
   bool accept = false;
   std::string scenario;
   
   statistics.set("funnel width", this->get_funnel_width());
   DEBUG << "\t\t" <<*this << '\n';

   bool funnel_acceptable = this->is_infeasibility_acceptable_to_funnel(trial_progress.infeasibility);

   if (funnel_acceptable)
   {
      if (this->in_restoration_phase)
      {
         if (this->armijo_sufficient_decrease(predicted_reduction.infeasibility, current_progress.infeasibility - trial_progress.infeasibility))
         {
            DEBUG << "Trial iterate (h-type) was accepted by satisfying the Armijo condition\n";
            accept = true;
         }
         else
         {
            DEBUG << "Trial iterate (h-type) was rejected by violating the Armijo condition\n";
         }
         scenario = "h-type Armijo";
         Iterate::number_eval_objective--;
      }
      else
      {
         // in filter and funnel methods, we construct an unconstrained measure by ignoring infeasibility and scaling the objective measure by 1
         const double current_merit = FunnelMethod::unconstrained_merit_function(current_progress);
         const double trial_merit = FunnelMethod::unconstrained_merit_function(trial_progress);
         const double merit_predicted_reduction = FunnelMethod::unconstrained_merit_function(predicted_reduction);
         DEBUG << "Current: (infeas., objective+auxiliary) = (" << current_progress.infeasibility << ", " << current_merit << ")\n";
         DEBUG << "Trial:   (infeas., objective+auxiliary) = (" << trial_progress.infeasibility << ", " << trial_merit << ")\n";
         DEBUG << "Unconstrained predicted reduction: " << merit_predicted_reduction << '\n';

         if (this->switching_condition(merit_predicted_reduction, current_progress.infeasibility, this->parameters.delta))
         {
            DEBUG << "\t\tTrial iterate satisfies switching condition ....\n";
            // unconstrained Armijo sufficient decrease condition (predicted reduction should be positive)
            const double objective_actual_reduction = this->compute_actual_objective_reduction(current_merit, current_progress.infeasibility,
                     trial_merit);
            DEBUG << "\t\tActual reduction: " << objective_actual_reduction << '\n';
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
            scenario = "current point";
         }
      }
   }
   else
   {
      DEBUG << "\t\tTrial iterate REJECTED by violating Funnel condition\n";
      scenario = "funnel";
   }
   statistics.set("status", std::string(accept ? "accepted" : "rejected") + " (" + scenario + ")");
   DEBUG << '\n';
   return accept;
}

// print the current funnel parameter
std::ostream& operator<<(std::ostream& stream, FunnelMethod& funnel) {
   stream << "************\n";
   stream << "\t\t  Current funnel width:\n";
   stream << "\t\t\t" << funnel.funnel_width << '\n';
   stream << "\t\t************\n";
   return stream;
}