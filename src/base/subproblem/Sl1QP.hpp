#ifndef Sl1QP_H
#define Sl1QP_H

#include <l1Relaxation.hpp>
#include "Subproblem.hpp"
#include "QPSolver.hpp"
#include "HessianEvaluation.hpp"
#include "ActiveSetMethod.hpp"
#include "l1Relaxation.hpp"

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
   Sl1QP(const Problem& problem, const std::string& QP_solver, const std::string& hessian_evaluation_method, bool use_trust_region, bool
   scale_residuals,
         double initial_parameter);

   void generate(const Problem& problem, const Iterate& current_iterate, double objective_multiplier, double trust_region_radius) override;
   void update_objective_multiplier(const Problem& problem, const Iterate& current_iterate, double objective_multiplier) override;

   Direction compute_direction(const Problem& problem, Iterate& current_iterate, double trust_region_radius) override;
   //Direction restore_feasibility(const Problem& problem, Iterate& current_iterate, Direction& phase_2_direction, double trust_region_radius)
   //override;
   int get_hessian_evaluation_count() override;

protected:
   /* use pointers to allow polymorphism */
   std::unique_ptr<QPSolver> solver; /*!< Subproblem solver */
   std::unique_ptr<HessianEvaluation> hessian_evaluation; /*!< Strategy to compute or approximate the Hessian */
   double penalty_parameter;
   Sl1QPParameters parameters;
   int number_variables;

   Sl1QP(const Problem& problem, const std::string& QP_solver, const std::string& hessian_evaluation_method, bool use_trust_region, bool
   scale_residuals, double initial_parameter, int number_variables);

   /* problem reformulation with elastic variables */
   // constraints l <= c(x) = u are reformulated as l <= c(x) - p + n <= u and p, n >= 0
   ElasticVariables elastic_variables_;

   size_t count_elastic_variables_(const Problem& problem);
   void set_variables_bounds_(const Problem& problem, const Iterate& current_iterate, double trust_region_radius) override;
   Direction solve_l1qp_subproblem_(const Problem& problem, Iterate& current_iterate, double trust_region_radius, double penalty_parameter);
   Direction compute_l1qp_step_(const Problem& problem, QPSolver& solver, Iterate& current_iterate, ConstraintPartition& constraint_partition,
         std::vector<double>& initial_solution, double trust_region_radius);
   Direction compute_l1qp_step_(const Problem& problem, QPSolver& solver, Iterate& current_iterate, double penalty_parameter,
         ElasticVariables& elastic_variables, double trust_region_radius);
   double compute_predicted_reduction_(const Problem& problem, Iterate& current_iterate, Direction& direction, double step_length);
   double compute_linearized_constraint_residual_(std::vector<double>& direction);
   double compute_error_(const Problem& problem, Iterate& iterate, Multipliers& multipliers, double penalty_parameter);
   double compute_complementarity_error_(const Problem& problem, Iterate& iterate, const Multipliers& multipliers) const override;
};

#endif // Sl1QP_H
