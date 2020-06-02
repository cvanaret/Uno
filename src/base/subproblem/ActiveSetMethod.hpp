#ifndef ACTIVESETMETHOD_H
#define ACTIVESETMETHOD_H

#include <memory>
#include "Subproblem.hpp"
#include "QPSolver.hpp"
#include "HessianEvaluation.hpp"

struct ElasticVariables {
    std::map<int, int> positive;
    std::map<int, int> negative;
};

class ActiveSetMethod : public Subproblem {
public:
    ActiveSetMethod(Problem& problem, bool scale_residuals);

    virtual Iterate evaluate_initial_point(Problem& problem, std::vector<double>& x, Multipliers& multipliers);
    
    void compute_optimality_measures(Problem& problem, Iterate& iterate);
    void compute_infeasibility_measures(Problem& problem, Iterate& iterate, SubproblemSolution& solution);

protected:
    /* QP subproblems */
    // optimality QP
    SubproblemSolution compute_qp_step_(Problem& problem, std::shared_ptr<QPSolver> solver, Iterate& current_iterate, double trust_region_radius);
    double compute_qp_predicted_reduction_(Iterate& current_iterate, SubproblemSolution& solution, double step_length);
    virtual std::vector<Range> generate_variables_bounds_(Problem& problem, Iterate& current_iterate, double trust_region_radius);
    // l1QP
    SubproblemSolution compute_l1qp_step_(Problem& problem, std::shared_ptr<QPSolver> solver, Iterate& current_iterate, ConstraintPartition& constraint_partition, std::vector<double>& initial_solution, double trust_region_radius);
    SubproblemSolution compute_l1qp_step_(Problem& problem, std::shared_ptr<QPSolver> solver, Iterate& current_iterate, double penalty_parameter, ElasticVariables& elastic_variables, double trust_region_radius);
    static void generate_elastic_variables_(Problem& problem, ElasticVariables& elastic_variables);
    void compute_l1_linear_objective_(Iterate& current_iterate, ConstraintPartition& constraint_partition);
    std::vector<double> generate_l1_multipliers_(Problem& problem, std::vector<double>& current_constraint_multipliers, ConstraintPartition& constraint_partition);
    std::vector<Range> generate_feasibility_bounds_(Problem& problem, std::vector<double>& current_constraints, ConstraintPartition& constraint_partition);
    static void recover_l1qp_active_set_(Problem& problem, SubproblemSolution& solution, const ElasticVariables& elastic_variables);
    
    /* LP subproblems */
    SubproblemSolution compute_lp_step_(Problem& problem, std::shared_ptr<QPSolver> solver, Iterate& current_iterate, double trust_region_radius);
    double compute_lp_predicted_reduction_(SubproblemSolution& solution, double step_length);
    SubproblemSolution compute_feasibility_lp_step_(Problem& problem, std::shared_ptr<QPSolver> solver, Iterate& current_iterate, SubproblemSolution& phase_II_solution, double trust_region_radius);
};

#endif // ACTIVESETMETHOD_H
