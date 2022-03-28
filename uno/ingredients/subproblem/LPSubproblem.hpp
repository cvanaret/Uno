#ifndef UNO_LPSUBPROBLEM_H
#define UNO_LPSUBPROBLEM_H

#include "ActiveSetSubproblem.hpp"
#include "solvers/QP/LPSolver.hpp"
#include "tools/Options.hpp"

class LPSubproblem : public ActiveSetSubproblem {
public:
   LPSubproblem(const Problem& problem, size_t max_number_variables, const Options& options);

   void build_objective_model(const Problem& problem, Iterate& current_iterate, double objective_multiplier) override;
   void build_constraint_model(const Problem& problem, Iterate& current_iterate) override;

   [[nodiscard]] Direction solve(Statistics& statistics, const Problem& problem, Iterate& current_iterate) override;
   [[nodiscard]] PredictedReductionModel generate_predicted_reduction_model(const Problem& problem, const Iterate& current_iterate,
         const Direction& direction) const override;
   [[nodiscard]] size_t get_hessian_evaluation_count() const override;
   [[nodiscard]] double get_proximal_coefficient() const override;

private:
   // use pointers to allow polymorphism
   const std::unique_ptr<LPSolver> solver; /*!< Solver that solves the subproblem */
};

#endif // UNO_LPSUBPROBLEM_H
