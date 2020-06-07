#include <cmath>
#include <map>
#include "SLPEQP.hpp"
#include "SLP.hpp"
#include "SQP.hpp"
#include "Constraint.hpp"
#include "Utils.hpp"
#include "Logger.hpp"
#include "QPSolverFactory.hpp"
#include "TrustRegion.hpp"
#include "LinearSolverFactory.hpp"

/* SLEQP: virtual class, implemented by SLEQP_TR and SLEQP_l2 */

SLPEQP::SLPEQP(Problem& problem, std::string LP_solver_name, std::string linear_solver_name, std::string hessian_evaluation_method, bool use_trust_region, bool scale_residuals):
ActiveSetMethod(problem, scale_residuals),
lp_solver(QPSolverFactory::create(LP_solver_name, problem.number_variables, problem.number_constraints, 0, false)),
linear_solver(LinearSolverFactory::create(linear_solver_name)),
hessian_evaluation(HessianEvaluationFactory::create(hessian_evaluation_method, problem.number_variables, false)) {
}

Direction SLPEQP::compute_step(Problem& problem, Iterate& current_iterate, double trust_region_radius) {
    /* compute first-order information */
    current_iterate.compute_objective_gradient(problem);
    current_iterate.compute_constraints(problem);
    current_iterate.compute_constraints_jacobian(problem);

    /***********/
    /* LP part */
    /***********/
    DEBUG << "SOLVING SLPEQP.LP\n";
    double LP_trust_region_radius = 10.;
    Direction direction_LP = this->compute_lp_step_(problem, this->lp_solver, current_iterate, LP_trust_region_radius);
    print_vector(DEBUG, direction_LP.x);
    
    /* set active trust region multipliers to 0 */
    TrustRegion::correct_active_set(direction_LP, LP_trust_region_radius);
    
    /************/
    /* EQP part */
    /************/
    DEBUG << "SOLVING SLPEQP.EQP\n";
    /* solve the EQP */
    if (direction_LP.status == INFEASIBLE) {
        return this->restore_feasibility(problem, current_iterate, direction_LP, trust_region_radius);
    }
    else {
        //SubproblemSolution solution_EQP = this->solver->solve_QP(problem.variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, current_iterate.hessian, solution_LP.x);
        //print_vector(DEBUG, solution_EQP.x);
        Direction direction_EQP = this->solve_eqp_(problem, current_iterate, direction_LP, trust_region_radius);
        direction_EQP.objective_multiplier = problem.objective_sign;
        direction_EQP.predicted_reduction = [&](double step_length) {
            return this->compute_qp_predicted_reduction_(current_iterate, direction_EQP, step_length);
        };
        return direction_EQP;
    }
}

Direction SLPEQP::solve_eqp_(Problem& problem, Iterate& current_iterate, Direction& phase_2_direction, double trust_region_radius) {
    // hs016 example
    DEBUG << "\nCurrent point: "; print_vector(DEBUG, current_iterate.x);
    // bounds
    std::vector<int> inactive_bound{1};
    std::set<int> active_bound_lb{0};
    std::set<int> active_bound_ub;
    
    // constraints
    std::vector<int> inactive_constraint{0};
    std::set<int> active_constraint_lb{1};
    std::set<int> active_constraint_ub;
    
    current_iterate.compute_constraints_jacobian(problem);
    current_iterate.compute_objective_gradient(problem);
    current_iterate.compute_constraints(problem);
    // use the multipliers from the LP solution (correct activity)
    this->hessian_evaluation->compute(problem, current_iterate, problem.objective_sign, phase_2_direction.multipliers.constraints);
    
    /* KKT matrix */
    
    int dimension = 2;
    COOMatrix kkt_matrix(dimension, 1);
    kkt_matrix.insert(current_iterate.hessian.matrix[2], 0, 0);
    kkt_matrix.insert(current_iterate.constraints_jacobian[1][1], 0, 1); // reduced Jacobian
    DEBUG << "KKT matrix:\n" << kkt_matrix;
    
    // TODO control inertia

    /* RHS */
    int rhs_size = inactive_bound.size() + active_constraint_lb.size() + active_constraint_ub.size();
    std::vector<double> rhs(rhs_size);
    // objective gradient: select inactive variables
    for (unsigned int i = 0; i < inactive_bound.size(); i++) {
        int original_i = inactive_bound[i];
        try {
            rhs[i] += -current_iterate.objective_gradient.at(original_i);
        }
        catch (const std::out_of_range& e) {}
    }
    // active constraints
    //TODO
    int active_constraint_index = 1;
    rhs[1] = -current_iterate.constraints[active_constraint_index] + problem.constraint_bounds[active_constraint_index].lb;
    DEBUG << "RHS: "; print_vector(DEBUG, rhs);
    
    this->linear_solver->factorize(kkt_matrix);
    this->linear_solver->solve(rhs);
    std::vector<double>& solution_EQP = rhs;
    DEBUG << "LP solution: "; print_vector(DEBUG, phase_2_direction.x);
    DEBUG << "Infeasible? " << (phase_2_direction.status == INFEASIBLE) << "\n";
    DEBUG << "EQP solution: "; print_vector(DEBUG, solution_EQP);
    
    Direction direction(phase_2_direction);
    direction.x[1] = solution_EQP[0];
    direction.multipliers.constraints[active_constraint_index] = -solution_EQP[1];
    return direction;
}

Direction SLPEQP::restore_feasibility(Problem& problem, Iterate& current_iterate, Direction& phase_2_direction, double trust_region_radius) {
    // TODO
    throw std::out_of_range("SLPEQP::restore_feasibility not implemented");
    return this->compute_qp_step_(problem, this->lp_solver, current_iterate, trust_region_radius);
}

///* SLPEQP_TR */
//
//SLPEQP_TR::SLPEQP_TR(Problem& problem, std::string LP_solver_name, std::string QP_solver_name, std::string hessian_evaluation_method, bool scale_residuals):
//SLPEQP(problem, LP_solver_name, hessian_evaluation_method, scale_residuals),
//qp_solver(QPSolverFactory::create(QP_solver_name, problem.number_variables, problem.number_constraints, problem.hessian_maximum_number_nonzeros + problem.number_variables, true)) {
//    this->hessian_evaluation->convexify = false;
//}
//
///* SLPEQP_l2 */
//
//SLPEQP_l2::SLPEQP_l2(Problem& problem, std::string LP_solver_name, std::string /*linear_solver_name*/, std::string hessian_evaluation_method, bool scale_residuals):
//SLPEQP(problem, LP_solver_name, hessian_evaluation_method, scale_residuals) {
//    this->hessian_evaluation->convexify = true;
//}
//
