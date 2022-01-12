#ifndef UNO_CONSTRAINTRELAXATIONSTRATEGY_H
#define UNO_CONSTRAINTRELAXATIONSTRATEGY_H

#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem/Direction.hpp"
#include "ingredients/subproblem/PredictedReductionModel.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "optimization/Problem.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Statistics.hpp"
#include "tools/Options.hpp"

struct ElasticVariables {
   SparseVector<size_t> positive;
   SparseVector<size_t> negative;
   explicit ElasticVariables(size_t capacity): positive(capacity), negative(capacity) {}
   [[nodiscard]] size_t size() const { return this->positive.size() + this->negative.size(); }
};

class ConstraintRelaxationStrategy {
public:
   ConstraintRelaxationStrategy(const Problem& problem, const Options& options);
   virtual ~ConstraintRelaxationStrategy() = default;

   virtual void initialize(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& first_iterate) = 0;
   virtual void create_current_subproblem(const Problem& problem, const Scaling& scaling, Iterate& current_iterate, double trust_region_radius) = 0;

   virtual Direction compute_feasible_direction(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& current_iterate) = 0;
   virtual Direction solve_feasibility_problem(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& current_iterate,
         const std::optional<std::vector<double>>& optional_phase_2_primal_direction,
         const std::optional<ConstraintPartition>& optional_constraint_partition) = 0;
   virtual Direction compute_second_order_correction(const Problem& problem, Iterate& trial_iterate);

   virtual bool is_acceptable(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, PredictedReductionModel& predicted_reduction_model, double step_length) = 0;
   virtual void register_accepted_iterate(Iterate& iterate);

   [[nodiscard]] PredictedReductionModel generate_predicted_reduction_model(const Problem& problem, const Direction& direction) const;
   [[nodiscard]] size_t get_hessian_evaluation_count() const;
   [[nodiscard]] size_t get_number_subproblems_solved() const;
   [[nodiscard]] SecondOrderCorrection soc_strategy() const;

protected:
   std::unique_ptr<Subproblem> subproblem;
   // possible problem reformulation with elastic variables. Constraints l <= c(x) <= u are reformulated as l <= c(x) - p + n <= u
   ElasticVariables elastic_variables;
   const double elastic_objective_coefficient;
   const size_t number_subproblem_variables;

   static size_t count_elastic_variables(const Problem& problem, bool subproblem_uses_slacks);
   static void generate_elastic_variables(const Problem& problem, ElasticVariables& elastic_variables, size_t number_variables,
         bool subproblem_uses_slacks);
   void add_elastic_variables_to_subproblem(Iterate& current_iterate);
   void remove_elastic_variables_from_subproblem();
   void remove_elastic_variables_from_direction(const Problem& problem, Direction& direction);
   void recover_active_set(const Problem& problem, Direction& direction);

public:
   const size_t max_number_subproblem_variables;
   const size_t number_constraints;
};

#endif //UNO_CONSTRAINTRELAXATIONSTRATEGY_H
