#ifndef UNO_ACTIVESETSUBPROBLEM_H
#define UNO_ACTIVESETSUBPROBLEM_H

#include "Subproblem.hpp"

class ActiveSetSubproblem : public Subproblem {
public:
   ActiveSetSubproblem(size_t max_number_variables, size_t number_constraints, SecondOrderCorrection soc_strategy, bool is_second_order_method,
         Norm residual_norm);
   ~ActiveSetSubproblem() override = default;

   void set_initial_point(const std::optional<std::vector<double>>& optional_initial_point) override;
   void set_elastic_variables(const l1RelaxedProblem& problem, Iterate& current_iterate) override;
   [[nodiscard]] double compute_optimality_measure(const NonlinearProblem& problem, Iterate& iterate) override;

protected:
   std::vector<double> initial_point;
   std::vector<Range> variable_displacement_bounds;
   std::vector<Range> linearized_constraint_bounds;

   void set_variable_displacement_bounds(const NonlinearProblem& problem, const Iterate& current_iterate);
   void set_linearized_constraint_bounds(const NonlinearProblem& problem, const std::vector<double>& current_constraints);
   static void compute_dual_displacements(const NonlinearProblem& problem, const Iterate& current_iterate, Direction& direction);
};

#endif // UNO_ACTIVESETSUBPROBLEM_H