#ifndef SLP_H
#define SLP_H

#include "Subproblem.hpp"
#include "solvers/QP/LPSolver.hpp"

class SLP : public Subproblem {
public:
   SLP(const Problem& problem, size_t max_number_variables, size_t number_constraints, const std::string& LP_solver_name);

   void create_current_subproblem(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) override;
   void build_objective_model(const Problem& problem, Iterate& current_iterate, double objective_multiplier) override;
   void set_initial_point(const std::vector<double>& point) override;
   Direction solve(Statistics& statistics, const Problem& problem, Iterate& current_iterate) override;
   [[nodiscard]] PredictedReductionModel generate_predicted_reduction_model(const Problem& problem, const Direction& direction) const override;
   [[nodiscard]] size_t get_hessian_evaluation_count() const override;

private:
   /* use pointers to allow polymorphism */
   const std::unique_ptr<LPSolver> solver; /*!< Solver that solves the subproblem */
   std::vector<double> initial_point;
};

#endif // SLP_H
