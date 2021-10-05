#ifndef GLOBALIZATIONMECHANISM_H
#define GLOBALIZATIONMECHANISM_H

#include "ingredients/constraint_relaxation/ConstraintRelaxationStrategy.hpp"
#include "optimization/Problem.hpp"
#include "tools/Statistics.hpp"

/*! \class GlobalizationMechanism
 * \brief Step control strategy
 *
 *  Strategy that promotes global convergence
 */
class GlobalizationMechanism {
public:
   GlobalizationMechanism(ConstraintRelaxationStrategy& constraint_relaxation_strategy, int max_iterations);
   virtual ~GlobalizationMechanism() = default;

   virtual void initialize(Statistics& statistics, const Problem& problem, Iterate& first_iterate) = 0;
   virtual std::tuple<Iterate, double> compute_acceptable_iterate(Statistics& statistics, const Problem& problem, Iterate& current_iterate) = 0;

   [[nodiscard]] int get_hessian_evaluation_count() const;
   [[nodiscard]] int get_number_subproblems_solved() const;

protected:
   /* references to allow polymorphism */
   ConstraintRelaxationStrategy& relaxation_strategy;
   int max_iterations; /*!< Maximum number of iterations */
   int number_iterations{0}; /*!< Current number of iterations */

   static Iterate assemble_trial_iterate(Iterate& current_iterate, Direction& direction, double step_length);
   static void print_warning(const char* message);
};

#endif // GLOBALIZATIONMECHANISM_H
