// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include "FunnelMethod.hpp"

FunnelMethod::FunnelMethod(bool in_restoration_phase, const Options& options) :
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
         options.get_double("funnel_beta"),
         options.get_double("funnel_gamma")
      }),
      in_restoration_phase(in_restoration_phase) {
}

void FunnelMethod::initialize(Statistics& /*statistics*/, const Iterate& initial_iterate, const Options& /*options*/) {
   // set the funnel upper bound
   double upper_bound = std::max(this->parameters.kappa_initial_upper_bound,
                                 this->parameters.kappa_initial_multiplication * initial_iterate.progress.infeasibility);

   // std::cout << "Initial kappa upper bound: " << upper_bound << std::endl;
   // std::cout << "Initial infeasibility: " << initial_iterate.progress.infeasibility << std::endl;
   // std::cout << "Initial funnel upper bound: " << upper_bound << std::endl;
   this->initial_funnel_upper_bound = upper_bound;
   this->funnel_width = this->initial_funnel_upper_bound;
//    this->current_iterate_acceptable_to_funnel = true;

}

void FunnelMethod::reset() {
   // re-initialize the restoration funnel
   // this->funnel->reset();
   // this->funnel->initial_upper_bound = this->initial_funnel_upper_bound;
}

void FunnelMethod::register_current_progress(const ProgressMeasures& /*current_progress_measures*/) {
   // const double current_optimality_measure = current_progress_measures.optimality(1.) + current_progress_measures.auxiliary_terms;
   // this->funnel->add(current_progress_measures.infeasibility, current_optimality_measure);
}

double FunnelMethod::get_infeasibility_upper_bound() const {
   return this->funnel_width;
}

void FunnelMethod::set_infeasibility_upper_bound(double new_upper_bound) {
   this->update_funnel_width(0.0, new_upper_bound);
}

bool FunnelMethod::is_infeasibility_acceptable_to_funnel(double infeasibility_measure) const {
   if (infeasibility_measure <= this->parameters.beta*this->funnel_width){
      return true;
   }
   else {
      DEBUG << "\t\tNot acceptable to funnel.\n";
      return false;
   }
}

bool FunnelMethod::is_infeasibility_acceptable(const ProgressMeasures& /*current_progress*/, const ProgressMeasures& /*trial_progress*/) const {
   // trial_iterate.evaluate_constraints(model);
   // const double trial_infeasibility = model.constraint_violation(trial_iterate.evaluations.constraints, progress_norm);
   // return this->is_infeasibility_acceptable_to_funnel(trial_infeasibility);
   return true;
}

bool FunnelMethod::switching_condition(double predicted_reduction, double current_infeasibility, double switching_fraction) const {
   return predicted_reduction > switching_fraction * std::pow(current_infeasibility, this->parameters.switching_infeasibility_exponent);
}

void FunnelMethod::update_funnel_width(double current_infeasibility_measure, double trial_infeasibility_measure) {

   this->funnel_width = std::max(this->parameters.kappa_infeasibility_1 *this->funnel_width, 
      trial_infeasibility_measure + this->parameters.kappa_infeasibility_2 * (current_infeasibility_measure - trial_infeasibility_measure));

   DEBUG << "\t\tNew funnel parameter is: " << this->funnel_width << "\n"; 
   
}

void FunnelMethod::update_funnel_width_restoration(double current_infeasibility_measure, double trial_infeasibility_measure) {

   if (trial_infeasibility_measure <= this->funnel_width){
      if (current_infeasibility_measure > this->funnel_width){
         this->funnel_width = std::min(this->parameters.kappa_infeasibility_1 *this->funnel_width,
         trial_infeasibility_measure + this->parameters.kappa_infeasibility_2 * (this->funnel_width - trial_infeasibility_measure));
         // std::cout << "Funnel Update Restoration: Current iterate outside funnel" << std::endl;
      } else {
         this->funnel_width = std::max(this->parameters.kappa_infeasibility_1 *this->funnel_width, 
         trial_infeasibility_measure + this->parameters.kappa_infeasibility_2 * (current_infeasibility_measure - trial_infeasibility_measure));
         // std::cout << "Funnel Update Restoration: Current iterate inside funnel" << std::endl;
      }
   } //else: do not reduce the funnel

   DEBUG << "\t\tNew funnel parameter is: " << this->funnel_width << "\n"; 
   
}

double FunnelMethod::get_funnel_width(){
   return this->funnel_width;
}

double FunnelMethod::compute_actual_reduction(double current_optimality_measure, double /*current_infeasibility_measure*/, double trial_optimality_measure) {
   return current_optimality_measure - trial_optimality_measure;
}

//! print: print the current funnel parameter
std::ostream& operator<<(std::ostream& stream, FunnelMethod& funnel) {
   stream << "************\n";
   stream << "\t\t  Current funnel width:\n";
   stream << "\t\t\t" << funnel.funnel_width << '\n';
   stream << "\t\t************\n";
   return stream;
}


/* check acceptability of step(s) (funnel & sufficient reduction)
 * funnel methods enforce an *unconstrained* sufficient decrease condition
 * precondition: feasible step
 * */
// bool FunnelMethod::is_iterate_acceptable(Statistics& /*statistics*/, const ProgressMeasures& /*current_progress*/,
//          const ProgressMeasures& /*trial_progress*/, const ProgressMeasures& /*predicted_reduction*/, double /*objective_multiplier*/) {
   // const double current_optimality_measure = current_progress_measures.optimality(1.) + current_progress_measures.auxiliary_terms;
   // const double trial_optimality_measure = trial_progress_measures.optimality(1.) + trial_progress_measures.auxiliary_terms;
   
   // // unconstrained predicted reduction:
   // // - ignore the predicted infeasibility reduction
   // // - scale the scaled optimality measure with 1
   // const double unconstrained_predicted_reduction = predicted_reduction.optimality(1.) + predicted_reduction.auxiliary_terms;
   // DEBUG << "\t\tCurrent: η = " << current_progress_measures.infeasibility << ",\t ω = " << current_optimality_measure << '\n';
   // DEBUG << "\t\tTrial:   η = " << trial_progress_measures.infeasibility << ",\t ω = " << trial_optimality_measure << '\n';
   // DEBUG << "\t\tUnconstrained predicted reduction: " << predicted_reduction.optimality(1.) << " + " << predicted_reduction.auxiliary_terms <<
   //       " = " <<  unconstrained_predicted_reduction << '\n';

   // GlobalizationStrategy::check_finiteness(current_progress_measures, 1.);
   // GlobalizationStrategy::check_finiteness(trial_progress_measures, 1.);
   // statistics.add_statistic("funnel width", this->funnel->get_funnel_size());
   
   // DEBUG << "\t\t" <<*this->funnel << '\n';

   // bool accept = false;
   // bool funnel_reduction_mechanism = false;
   // bool funnel_acceptable = false;

   // funnel_acceptable = this->is_infeasibility_acceptable(trial_progress_measures.infeasibility);
   
   // // check acceptance   
   // if (funnel_acceptable) {
   //    DEBUG << "\t\tFunnel condition acceptable or in feasibility restoration phase\n";

   //       const double actual_reduction = this->funnel->compute_actual_reduction(current_optimality_measure, current_progress_measures.infeasibility,
   //             trial_optimality_measure);
   //       DEBUG << "\t\tActual reduction: " << actual_reduction << '\n';

   //       DEBUG << "\t\tCheck the switching condition ....\n";
   //       // switching condition: the unconstrained predicted reduction is sufficiently positive
   //       if (this->switching_condition(unconstrained_predicted_reduction, current_progress_measures.infeasibility, this->parameters.delta)) {
   //          // unconstrained Armijo sufficient decrease condition (predicted reduction should be positive)
   //          DEBUG << "\t\tCheck the armijo condition for descent ....\n";
   //          if (this->armijo_sufficient_decrease(unconstrained_predicted_reduction, actual_reduction)) {
   //             DEBUG << "\t\tTrial iterate was ACCEPTED by satisfying Armijo condition\n";
   //             accept = true;
   //             // decrease funnel here ......
   //          }
   //          else { // switching condition holds, but not Armijo condition
   //             DEBUG << "\t\tArmijo condition not satisfied, trial iterate REJECTED\n";
   //          }
   //       }
   //       else {
   //          DEBUG << "\t\tTrial iterate violates switching condition ...\n";
   //          funnel_reduction_mechanism = true;
   //          accept = true; // accept the step and reduce the tr-radius
   //       }

   // } else {
   //    funnel_reduction_mechanism = false;
   //    // step is rejected
   //    DEBUG << "\t\tFunnel condition NOT acceptable\n";
   // }

   // if (funnel_reduction_mechanism){
   //     DEBUG << "\t\tEntering funnel reduction mechanism\n";

   //    // // Feasibility measures
   //    const double current_infeasibility_measure = current_progress_measures.infeasibility;
   //    const double trial_infeasibility_measure = trial_progress_measures.infeasibility;
   //    this->funnel->update_funnel_parameter(current_infeasibility_measure, 
   //                                           trial_infeasibility_measure);
   //    this->funnel_width = this->funnel->get_funnel_size()                                             ;
   // }

   // if (accept){
   //    if (trial_progress_measures.infeasibility < this->funnel->get_funnel_size()){
   //       this->current_iterate_acceptable_to_funnel = true;
   //    }
   //    if (this->current_phase == 1){
   //       //update funnel ....
   //    }
   // }

   // return accept;
//    return true;
// }

bool FunnelMethod::is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
         const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction, double /*objective_multiplier*/) {
   const double current_optimality_measure = current_progress.objective(1.) + current_progress.auxiliary;
   const double trial_optimality_measure = trial_progress.objective(1.) + trial_progress.auxiliary;
   
   const double current_infeasibility_measure = current_progress.infeasibility;
   const double trial_infeasibility_measure = trial_progress.infeasibility;

   // unconstrained predicted reduction:
   // - ignore the predicted infeasibility reduction
   // - scale the scaled optimality measure with 1
   const double unconstrained_predicted_reduction = predicted_reduction.objective(1.) + predicted_reduction.auxiliary;

   DEBUG << "\t\tCurrent: η = " << current_progress.infeasibility << ",\t ω = " << current_optimality_measure << '\n';
   DEBUG << "\t\tTrial:   η = " << trial_progress.infeasibility << ",\t ω = " << trial_optimality_measure << '\n';
   DEBUG << "\t\tUnconstrained predicted reduction: " << predicted_reduction.objective(1.) << " + " << predicted_reduction.auxiliary <<
         " = " <<  unconstrained_predicted_reduction << '\n';
   DEBUG << "\t\tUnconstrained predicted infeasibility reduction: " << predicted_reduction.infeasibility << '\n';

   statistics.set("funnel width", this->get_funnel_width());
   
   DEBUG << "\t\t" <<*this << '\n';

   bool accept = false;
   bool funnel_reduction_mechanism = false;
   bool funnel_acceptable = this->is_infeasibility_acceptable_to_funnel(trial_progress.infeasibility);

   if (funnel_acceptable) {
      if (this->switching_condition(unconstrained_predicted_reduction, current_progress.infeasibility, this->parameters.delta)) {
         DEBUG << "\t\tTrial iterate satisfies switching condition ....\n";
         // check acceptance 
         const double actual_reduction = this->compute_actual_reduction(current_optimality_measure, current_progress.infeasibility,
               trial_optimality_measure);
         DEBUG << "\t\tActual reduction: " << actual_reduction << '\n';
         
         // unconstrained Armijo sufficient decrease condition (predicted reduction should be positive)
         if (this->armijo_sufficient_decrease(unconstrained_predicted_reduction, actual_reduction)) {
            DEBUG << "\t\tTrial iterate was ACCEPTED by satisfying Armijo condition\n";
            accept = true;
            statistics.set("status", "accepted (Armijo O)");
            if (this->in_restoration_phase){
               DEBUG << "\t\tEntering funnel reduction mechanism in Restoration Phase\n";
               this->update_funnel_width_restoration(current_optimality_measure, trial_optimality_measure);
            }
         }
         else { // switching condition holds, but not Armijo condition
            DEBUG << "\t\tArmijo condition not satisfied, trial iterate REJECTED\n";
            statistics.set("status", "rejected (Armijo)");
         }
      } else {
         DEBUG << "\t\tTrial iterate ACCEPTED by violating the switching condition ...\n";
         funnel_reduction_mechanism = true;
         accept = true; // accept the step and reduce the tr-radius
         statistics.set("status", "accepted (!switching)");
      }
   } else {
      DEBUG << "\t\tTrial iterate REJECTED by violating Funnel condition\n";
      statistics.set("status", "rejected (funnel)");
      // acceptable = false;
      // funnel_reduction_mechanism = false;
   }

   if (funnel_reduction_mechanism){ // steps needs to be accepted for this...
       DEBUG << "\t\tEntering funnel reduction mechanism\n";

      // // Feasibility measures
      this->update_funnel_width(current_infeasibility_measure, 
                                             trial_infeasibility_measure);
      // this->funnel_width = this->get_funnel_width(); //?                                            ;
   }

   return accept;
}
