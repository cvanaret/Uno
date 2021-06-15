#ifndef GLOBALIZATIONMECHANISM_H
#define GLOBALIZATIONMECHANISM_H

#include "Problem.hpp"
#include "GlobalizationStrategy.hpp"
#include "Statistics.hpp"

/*! \class GlobalizationMechanism
 * \brief Step control strategy
 *
 *  Strategy that promotes global convergence
 */
class GlobalizationMechanism {
public:
   GlobalizationMechanism(ConstraintRelaxationStrategy& constraint_relaxation_strategy, int max_iterations);
   virtual ~GlobalizationMechanism() = default;

   virtual Iterate initialize(Statistics& statistics, const Problem& problem, std::vector<double>& x, Multipliers& multipliers) = 0;
   virtual std::pair<Iterate, Direction> compute_acceptable_iterate(Statistics& statistics, const Problem& problem, Iterate& current_iterate) = 0;

   [[nodiscard]] int get_hessian_evaluation_count() const;
   [[nodiscard]] int get_number_subproblems_solved() const;

protected:
   /* references to allow polymorphism */
   ConstraintRelaxationStrategy& relaxation_strategy;
   int max_iterations; /*!< Maximum number of iterations */
   int number_iterations; /*!< Current number of iterations */
   // preallocated vector to receive the trial primal variables
   std::vector<double> trial_primals_;

   Iterate assemble_trial_iterate(const Iterate& current_iterate, const Direction& direction, double step_length);
   virtual void print_acceptance_(const Iterate& iterate);
   static void print_warning_(const char* message);
};

#endif // GLOBALIZATIONMECHANISM_H
