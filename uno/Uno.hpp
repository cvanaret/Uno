// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_H
#define UNO_H

#include "optimization/Model.hpp"
#include "optimization/TerminationStatus.hpp"
#include "ingredients/globalization_mechanism/GlobalizationMechanism.hpp"

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
   Uno(GlobalizationMechanism& globalization_mechanism, const Options& options);

   [[nodiscard]] Result solve(const Model& model, Iterate& first_iterate, const Options& options);

private:
   GlobalizationMechanism& globalization_mechanism; /*!< Step control strategy (trust region or line-search) */
   const double tolerance; /*!< Tolerance of the termination criteria */
   const size_t max_iterations; /*!< Maximum number of iterations */
   const double small_step_factor{100.};

   static Statistics create_statistics(const Model& model, const Options& options);
   static void add_statistics(Statistics& statistics, const Model& model, const Iterate& new_iterate, size_t major_iterations);
   [[nodiscard]] bool termination_criterion(TerminationStatus current_status, size_t iteration) const;
   [[nodiscard]] TerminationStatus check_termination(const Model& model, Iterate& current_iterate, double step_norm) const;
};

#endif // UNO_H
