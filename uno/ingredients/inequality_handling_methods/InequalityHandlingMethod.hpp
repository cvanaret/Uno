// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INEQUALITYHANDLINGMETHOD_H
#define UNO_INEQUALITYHANDLINGMETHOD_H

#include <memory>
#include <string>
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/regularization_strategies/RegularizationStrategy.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   // forward declarations
   class Direction;
   class Iterate;
   class l1RelaxedProblem;
   class Model;
   class Multipliers;
   class OptimizationProblem;
   class Options;
   class Statistics;
   template <typename IndexType, typename ElementType>
   class SymmetricMatrix;
   template <typename ElementType>
   class Vector;
   class WarmstartInformation;
   
   class InequalityHandlingMethod {
   public:
      InequalityHandlingMethod(const std::string& hessian_model, const std::string& regularization_strategy, size_t dimension,
            size_t number_hessian_nonzeros, bool convexify, const Options& options);
      virtual ~InequalityHandlingMethod() = default;

      // virtual methods implemented by subclasses
      virtual void initialize_statistics(Statistics& statistics, const Options& options) = 0;
      virtual void generate_initial_iterate(const OptimizationProblem& problem, Iterate& initial_iterate) = 0;
      virtual void solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate, const Multipliers& current_multipliers,
            Direction& direction, WarmstartInformation& warmstart_information) = 0;

      void set_trust_region_radius(double new_trust_region_radius);
      virtual void initialize_feasibility_problem(const l1RelaxedProblem& problem, Iterate& current_iterate) = 0;
      virtual void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate) = 0;
      [[nodiscard]] virtual double proximal_coefficient(const Iterate& current_iterate) const = 0;
      virtual void exit_feasibility_problem(const OptimizationProblem& problem, Iterate& trial_iterate) = 0;

      // progress measures
      [[nodiscard]] virtual double hessian_quadratic_product(const Vector<double>& primal_direction) const = 0;
      virtual void set_auxiliary_measure(const Model& model, Iterate& iterate) = 0;
      [[nodiscard]] virtual double compute_predicted_auxiliary_reduction_model(const Model& model, const Iterate& current_iterate,
            const Vector<double>& primal_direction, double step_length) const = 0;

      virtual void postprocess_iterate(const OptimizationProblem& problem, Iterate& iterate) = 0;

      [[nodiscard]] size_t get_hessian_evaluation_count() const;
      virtual void set_initial_point(const Vector<double>& initial_point) = 0;

      size_t number_subproblems_solved{0};
      // when the parameterization of the subproblem (e.g. penalty or barrier parameter) is updated, signal it
      bool subproblem_definition_changed{false};

   protected:
      const std::unique_ptr<HessianModel> hessian_model; /*!< Strategy to evaluate or approximate the Hessian */
      const std::unique_ptr<RegularizationStrategy<double>> regularization_strategy;
      double trust_region_radius{INF<double>};
   };
} // namespace

#endif // UNO_INEQUALITYHANDLINGMETHOD_H