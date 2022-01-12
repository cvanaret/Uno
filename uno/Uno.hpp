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
   Scaling scaling;
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
    */
   Uno(GlobalizationMechanism& globalization_mechanism, const Options& options);

   /*!
    *  Solve a given problem with initial primal and dual variables
    */
   Result solve(const Problem& problem, Iterate& first_iterate, bool scale_functions, bool enforce_linear_constraints);

private:
   GlobalizationMechanism& globalization_mechanism; /*!< Step control strategy (trust region or line-search) */
   const double tolerance; /*!< Tolerance of the termination criteria */
   const size_t max_iterations; /*!< Maximum number of iterations */
   const double small_step_factor{100.};

   static Statistics create_statistics(const Problem& problem);
   static void add_statistics(Statistics& statistics, const Problem& problem, const Iterate& new_iterate, size_t major_iterations);
   [[nodiscard]] bool termination_criterion(TerminationStatus current_status, size_t iteration) const;
   [[nodiscard]] TerminationStatus check_termination(const Problem& problem, const Iterate& current_iterate, double step_norm) const;
   static void postsolve_solution(const Problem& problem, const Scaling& scaling, Iterate& current_iterate, TerminationStatus termination_status);
};

#endif // UNO_H
