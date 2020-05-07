#include "TubeStrategy.hpp"
#include "Utils.hpp"

/* TODO initial value of upper_bound? */
TubeStrategy::TubeStrategy(StepConstants& constants, Tolerances& tolerances):
	Strategy(constants, tolerances),
	tube_restoration(Tube(0., constants)),
	tube_optimality(Tube(0., constants)) {
}

bool TubeStrategy::check_step(Problem& problem, Iterate& current_iterate, Step& step, Phase& phase) {
	//int errflag = 0;
	//int iprint = 0;
	bool accept = false;

	///* infeasible QP => go from phase II (optimality) to I (restoration) */
	//if (phase == OPTIMALITY && (step.status == UNBOUNDED || step.status == INFEASIBLE)) {
		//phase = RESTORATION;
		///* add [h,f] (c/s violation) to filter, entering restoration */
		//// what should be the value of trial_iterate.residual?
		//// TODO this->tube_restoration.update(current_iterate.residual, trial_iterate.residual);
		
		///* compute values for current h=feasible constraints, f=infeasible constraints */
		//double residual = this->problem->feasible_residual_norm(step.feasible_constraints, step.infeasible_constraints, current_iterate.constraints);
		//double objective = this->problem->infeasible_residual_norm(step.infeasible_constraints, current_iterate.constraints);
		
		//this->tube_restoration.upper_bound = max(this->tube_optimality.upper_bound, residual + objective);
	//} 

	///* generate the trial point */
	//Iterate trial_iterate = this->generate_trial_iterate(current_iterate.x, step.d);
	
	//if (phase == OPTIMALITY || (phase == RESTORATION && step.status <= 1)) {
		//trial_iterate.objective = this->problem->objective.eval(trial_iterate.x);
	//}

	///* check whether we can change from phase I (restoration) to II (optimality) */
	//if (phase == RESTORATION && step.status <= 1) {
		///* check whether (trial_eval.residual, trial_eval.objective) acceptable to Tube */
		//bool acceptable = this->tube_optimality.query(trial_iterate.residual);
		///* QP feasible and (trial_eval.residual, trial_eval.objective) acceptable then change back to phase II */
		//if (acceptable) {
			//phase = OPTIMALITY;
			///* recompute [h,f] for phase 2 (optimization) */
			//current_iterate.objective = this->problem->objective.eval(current_iterate.x);
			
			//current_iterate.residual = this->problem->l1_inf_norm(current_iterate.constraints);
			
			//if (2 <= iprint) { 
				//fprintf(stdout,"   QP feasible & (h,f) acceptable -> return to phase II\n");
				//// TODO PrintExitRest(current_iterate.residual, *f);
			//}
		//}
		//else if (2 <= iprint) {
			//fprintf(stdout,"   QP feasible BUT (h,f) NOT acceptable -> remain in phase I\n");
		//}
	//}

	///* distinguish phase 1 (restoration) & phase 2 (optimization) */
	//if (phase == RESTORATION) {
		///* compute new objective value & constraint norm for RESTORATION */
		//residual = this->problem->feasible_residual_norm(step.feasible_constraints, step.infeasible_constraints, current_iterate.constraints);
		//objective = this->problem->infeasible_residual_norm(step.infeasible_constraints, current_iterate.constraints);
		///* check acceptability criteria */
		//accept = this->tube_restoration.query(trial_iterate.residual);
	//}
	//else {
		//accept = this->tube_optimality.query(trial_iterate.residual);
	//}

	//if (accept) {
		///* check sufficient reduction condition */
		//double actual_reduction = current_iterate.objective - trial_iterate.objective;
		//double square_residual = current_iterate.residual*current_iterate.residual;
		//if (actual_reduction >= this->constants.Sigma*step.predicted_reduction) {          // successful f-type step
			//accept = true;
		//}
		//else if (step.predicted_reduction > this->constants.Delta*square_residual) {       // unsuccessful h-type step
			//accept = false;
		//}
		///* add only successful h-type steps to tube */
		//if (accept && step.predicted_reduction < this->constants.Delta*square_residual) {
			//if (phase == RESTORATION) {
				//this->tube_restoration.update(current_iterate.residual, trial_iterate.residual);
			//}
			//else {
				//this->tube_optimality.update(current_iterate.residual, trial_iterate.residual);
			//}
		//}
	//}

	return accept;
}
