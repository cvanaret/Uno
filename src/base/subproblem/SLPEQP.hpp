#ifndef SLPEQP_H
#define SLPEQP_H

#include "Subproblem.hpp"
#include "QPSolver.hpp"
#include "MA57Solver.hpp"
#include "SQP.hpp"
#include "SLP.hpp"

/*! \class QPApproximation
 * \brief QP local approximation
 *
 *  Quadratic approximation
 */
class SLPEQP : public Subproblem {
public:
    /*!
     *  Constructor
     * 
     * \param solver: solver that solves the subproblem
     */
    SLPEQP(Problem& problem, std::string QP_solver_name, std::string hessian_evaluation_method, bool use_trust_region, bool scale_residuals);

    SubproblemSolution compute_step(Problem& problem, Iterate& current_iterate, double trust_region_radius=INFINITY) override;
    SubproblemSolution restore_feasibility(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution, double trust_region_radius=INFINITY) override;
    
    void compute_optimality_measures(Problem& problem, Iterate& iterate) override;
    void compute_infeasibility_measures(Problem& problem, Iterate& iterate, SubproblemSolution& solution) override;
    
    bool phase_1_required(SubproblemSolution& solution) override;

    /* use a reference to allow polymorphism */
    SLP lp_subproblem;
    SQP eqp_subproblem;

private:
    void fix_active_constraints(Problem& problem, ActiveSet& active_set, std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds);
    double compute_predicted_reduction(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length);
};

//class SLPEQP_l2 : public SLPEQP {
//public:
//    /*!
//     *  Constructor
//     * 
//     * \param solver: solver that solves the subproblem
//     */
//    SLPEQP_l2();
//
//    Iterate initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, int number_variables, int number_constraints, bool use_trust_region) override;
//
//    SubproblemSolution compute_optimality_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds) override;
//    SubproblemSolution compute_infeasibility_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds, SubproblemSolution& phase_II_solution) override;
//    void compute_measures(Problem& problem, Iterate& iterate) override;
//    bool phase_1_required(SubproblemSolution& solution) override;
//    double compute_predicted_reduction(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length) override;
//
//    /* use a reference to allow polymorphism */
//    MA57Solver solver; /*!< Solver that solves the subproblem */
//
//private:
//    void fix_active_constraints(Problem& problem, ActiveSet& active_set, std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds);
//};
//
//class SLPEQP_TR : public SLPEQP {
//public:
//    /*!
//     *  Constructor
//     * 
//     * \param solver: solver that solves the subproblem
//     */
//    SLPEQP_TR(QPSolver& solver);
//
//    Iterate initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, int number_variables, int number_constraints, bool use_trust_region) override;
//
//    SubproblemSolution compute_optimality_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds) override;
//    SubproblemSolution compute_infeasibility_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds, SubproblemSolution& phase_II_solution) override;
//    void compute_measures(Problem& problem, Iterate& iterate) override;
//    bool phase_1_required(SubproblemSolution& solution) override;
//    double compute_predicted_reduction(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length) override;
//
//    /* use a reference to allow polymorphism */
//    QPSolver& solver; /*!< Solver that solves the subproblem */
//
//private:
//    void fix_active_constraints(Problem& problem, ActiveSet& active_set, std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds);
//};

#endif // SLPEQP_H
