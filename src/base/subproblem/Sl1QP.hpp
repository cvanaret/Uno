#ifndef Sl1QP_H
#define Sl1QP_H

#include <Relaxation.hpp>
#include "Subproblem.hpp"
#include "QPSolver.hpp"
#include "HessianEvaluation.hpp"
#include "ActiveSetMethod.hpp"
#include "Relaxation.hpp"

struct Sl1QPParameters {
   double tau;
   double epsilon1;
   double epsilon2;
};

/*! \class QPApproximation
 * \brief QP local approximation
 *
 *  Quadratic approximation
 */
class Sl1QP : public ActiveSetMethod {
public:
   /*!
    *  Constructor
    *
    * \param solver: solver that solves the subproblem
    */
   Sl1QP(Problem& problem, std::string QP_solver, std::string hessian_evaluation_method, bool use_trust_region, bool scale_residuals,
         double initial_parameter);

   std::vector<Direction> compute_directions(Problem& problem, Iterate& current_iterate, double objective_multiplier,
         double trust_region_radius = INFINITY) override;
   std::vector<Direction> restore_feasibility(Problem& problem, Iterate& current_iterate, Direction& phase_2_direction,
         double trust_region_radius = INFINITY) override;

   /* use pointers to allow polymorphism */
   std::unique_ptr<QPSolver> solver; /*!< Subproblem solver */
   std::unique_ptr<HessianEvaluation> hessian_evaluation; /*!< Strategy to compute or approximate the Hessian */
   double penalty_parameter;
   Sl1QPParameters parameters;
   int number_variables;

protected:
   Sl1QP(Problem& problem, std::string QP_solver, std::string hessian_evaluation_method, bool use_trust_region, bool scale_residuals,
         double initial_parameter, int number_variables);

   /* problem reformulation with elastic variables */
   // constraints l <= c(x) = u are reformulated as c(x) - p + n
   ElasticVariables elastic_variables_;

   size_t count_elastic_variables_(Problem& problem);
   void evaluate_optimality_iterate_(Problem& problem, Iterate& current_iterate, double penalty_parameter);
   std::vector<Range> generate_variables_bounds_(Problem& problem, Iterate& current_iterate, double trust_region_radius) override;
   Direction solve_l1qp_subproblem_(Problem& problem, Iterate& current_iterate, double trust_region_radius, double penalty_parameter);
   double compute_predicted_reduction_(Problem& problem, Iterate& current_iterate, Direction& direction, double step_length);
   double compute_linearized_constraint_residual_(std::vector<double>& direction);
   double compute_error_(Problem& problem, Iterate& iterate, Multipliers& multipliers, double penalty_parameter);
   double compute_complementarity_error_(const Problem& problem, Iterate& iterate, const Multipliers& multipliers) const override;
};

#endif // Sl1QP_H
