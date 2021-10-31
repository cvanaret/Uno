#ifndef UNO_H
#define UNO_H

#include "optimization/Problem.hpp"
#include "ingredients/mechanism/GlobalizationMechanism.hpp"

enum TerminationStatus {
   NOT_OPTIMAL = 0,
   KKT_POINT, /* feasible stationary point */
   FJ_POINT, /* infeasible stationary point */
   FEASIBLE_SMALL_STEP,
   INFEASIBLE_SMALL_STEP
};

struct Result {
   Result() = delete;

   TerminationStatus status;
   Iterate solution;
   size_t number_variables;
   size_t number_constraints;
   size_t iteration;
   double cpu_time;
   size_t objective_evaluations;
   size_t constraint_evaluations;
   size_t jacobian_evaluations;
   size_t hessian_evaluations;
   size_t number_subproblems_solved;

   void print(bool print_solution) const;
};

/*! \class Uno
 * \brief Uno
 *
 *  UNO solver
 */
class Uno {
public:
   /*!
    *  Constructor
    *
    * \param globalization_strategy: strategy to promote global convergence
    * \param tolerance: tolerance for termination criteria
    */
   Uno(GlobalizationMechanism& globalization_mechanism, double tolerance, size_t max_iterations);

   /*!
    *  Solve a given problem with initial primal and dual variables
    *
    * \param problem: optimization problem
    * \param x: primal variables
    * \param multipliers: Lagrange multipliers/dual variables
    */
   Result solve(const Problem& problem, Iterate& first_iterate, bool use_preprocessing);

private:
   GlobalizationMechanism& globalization_mechanism; /*!< Step control strategy (trust region or line-search) */
   double tolerance; /*!< Tolerance of the termination criteria */
   size_t max_iterations; /*!< Maximum number of iterations */

   static Statistics create_statistics();
   static void add_statistics(Statistics& statistics, const Iterate& new_iterate, size_t major_iterations);
   [[nodiscard]] bool termination_criterion(TerminationStatus current_status, size_t iteration) const;
   TerminationStatus check_termination(const Problem& problem, Iterate& current_iterate, double step_norm) const;
};

#endif // UNO_H
