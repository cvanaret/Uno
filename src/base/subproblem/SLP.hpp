#ifndef SLP_H
#define SLP_H

#include "ActiveSetMethod.hpp"

class SLP : public ActiveSetMethod {
public:
    SLP(Problem& problem, std::string QP_solver_name, bool use_trust_region, bool scale_residuals);

    SubproblemSolution compute_step(Problem& problem, Iterate& current_iterate, double trust_region_radius = INFINITY) override;
    SubproblemSolution restore_feasibility(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution, double trust_region_radius = INFINITY) override;
    
private:
    void evaluate_optimality_iterate_(Problem& problem, Iterate& current_iterate);
    void evaluate_feasibility_iterate_(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution);
};

#endif // SLP_H
