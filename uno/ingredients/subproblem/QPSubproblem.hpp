#ifndef UNO_QPSUBPROBLEM_H
#define UNO_QPSUBPROBLEM_H

#include "ActiveSetSubproblem.hpp"
#include "HessianModel.hpp"
#include "ingredients/strategy/NonlinearProblem.hpp"
#include "solvers/QP/QPSolver.hpp"
#include "tools/Options.hpp"

class QPSubproblem : public ActiveSetSubproblem {
public:
   QPSubproblem(const NonlinearProblem& problem, size_t max_number_variables, const Options& options);

   [[nodiscard]] Direction solve(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate) override;
   [[nodiscard]] PredictedReductionModel generate_predicted_reduction_model(const NonlinearProblem& problem, const Direction& direction) const override;
   [[nodiscard]] size_t get_hessian_evaluation_count() const override;
   [[nodiscard]] double get_proximal_coefficient() const override;

protected:
   // use pointers to allow polymorphism
   const std::unique_ptr<QPSolver> solver; /*!< Solver that solves the subproblem */
   const double proximal_coefficient;

   // evaluations
   const std::unique_ptr<HessianModel> hessian_model; /*!< Strategy to evaluate or approximate the Hessian */
   SparseVector<double> objective_gradient;
   std::vector<double> constraints;
   std::vector<SparseVector<double>> constraint_jacobian;

   void evaluate_problem(const NonlinearProblem& problem, Iterate& current_iterate);
};

#endif // UNO_QPSUBPROBLEM_H
