// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_GLOBALIZATIONMECHANISM_H
#define UNO_GLOBALIZATIONMECHANISM_H

#include "ingredients/subproblem/Direction.hpp"

// forward declarations
class ConstraintRelaxationStrategy;
class Iterate;
class Model;
class Options;
class Statistics;

class GlobalizationMechanism {
public:
   explicit GlobalizationMechanism(ConstraintRelaxationStrategy& constraint_relaxation_strategy);
   virtual ~GlobalizationMechanism() = default;

   virtual void initialize(Statistics& statistics, Iterate& initial_iterate, const Options& options) = 0;
   virtual void compute_next_iterate(Statistics& statistics, const Model& model, Iterate& current_iterate, Iterate& trial_iterate) = 0;

   [[nodiscard]] size_t get_hessian_evaluation_count() const;
   [[nodiscard]] size_t get_number_subproblems_solved() const;

protected:
   // reference to allow polymorphism
   ConstraintRelaxationStrategy& constraint_relaxation_strategy; /*!< Constraint relaxation strategy */
   Direction direction;

   static void assemble_trial_iterate(const Model& model, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double primal_step_length, double dual_step_length);
};

#endif // UNO_GLOBALIZATIONMECHANISM_H
