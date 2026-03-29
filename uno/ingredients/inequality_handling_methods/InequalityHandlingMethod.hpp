// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INEQUALITYHANDLINGMETHOD_H
#define UNO_INEQUALITYHANDLINGMETHOD_H

#include <memory>
#include <string>

namespace uno {
   // forward declarations
   class Evaluations;
   class Iterate;
   class l1RelaxedProblem;
   class OptimizationProblem;
   class Parameterization;
   class Statistics;
   
   class InequalityHandlingMethod {
   public:
      InequalityHandlingMethod() = default;
      virtual ~InequalityHandlingMethod() = default;

      virtual void initialize_statistics(Statistics& statistics) = 0;
      [[nodiscard]] virtual std::unique_ptr<OptimizationProblem> reformulate(const OptimizationProblem& problem,
         Parameterization& parameterization) = 0;
      [[nodiscard]] virtual bool update_parameterization(Statistics& statistics, const OptimizationProblem& problem,
         const Iterate& current_iterate, Parameterization& parameterization) = 0;

      virtual void initialize_feasibility_problem(Iterate& current_iterate) = 0;
      virtual void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate, Evaluations& evaluations) = 0;
      [[nodiscard]] virtual double proximal_coefficient() const = 0;

      [[nodiscard]] virtual std::string get_name() const = 0;
   };
} // namespace

#endif // UNO_INEQUALITYHANDLINGMETHOD_H