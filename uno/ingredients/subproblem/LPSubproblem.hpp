#ifndef UNO_LPSUBPROBLEM_H
#define UNO_LPSUBPROBLEM_H

#include "Subproblem.hpp"
#include "solvers/QP/LPSolver.hpp"
#include "tools/Options.hpp"

class LPSubproblem : public Subproblem {
public:
   LPSubproblem(const Problem& problem, size_t max_number_variables, const Options& options);

   void build_current_subproblem(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) override;
   void build_objective_model(const Problem& problem, Iterate& current_iterate, double objective_multiplier) override;
   void set_initial_point(const std::vector<double>& point) override;
   [[nodiscard]] Direction solve(Statistics& statistics, const Problem& problem, Iterate& current_iterate) override;
   [[nodiscard]] PredictedReductionModel generate_predicted_reduction_model(const Problem& problem, const Direction& direction) const override;
   [[nodiscard]] size_t get_hessian_evaluation_count() const override;

private:
   // use pointers to allow polymorphism
   const std::unique_ptr<LPSolver> solver; /*!< Solver that solves the subproblem */
   std::vector<double> initial_point;
};

#endif // UNO_LPSUBPROBLEM_H
