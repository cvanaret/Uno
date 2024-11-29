// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "QPSubproblem.hpp"
#include "ingredients/local_models/LagrangeNewtonSubproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "preprocessing/Preprocessing.hpp"
#include "reformulation/OptimizationProblem.hpp"
#include "solvers/QPSolver.hpp"
#include "solvers/QPSolverFactory.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   QPSubproblem::QPSubproblem(size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros,
         size_t number_jacobian_nonzeros, size_t number_hessian_nonzeros, const Options& options) :
         InequalityConstrainedMethod(options.get_string("hessian_model"), number_variables, number_hessian_nonzeros,
               options.get_string("globalization_mechanism") == "LS", options),
         use_regularization(options.get_string("globalization_mechanism") != "TR" || options.get_bool("convexify_QP")),
         enforce_linear_constraints_at_initial_iterate(options.get_bool("enforce_linear_constraints")),
         // maximum number of Hessian nonzeros = number nonzeros + possible diagonal inertia correction
         solver(QPSolverFactory::create(number_variables, number_constraints, number_objective_gradient_nonzeros, number_jacobian_nonzeros,
               // if the QP solver is used during preprocessing, we need to allocate the Hessian with at least number_variables elements
               std::max(this->enforce_linear_constraints_at_initial_iterate ? number_variables : 0, hessian_model->hessian.capacity()),
               options)) {
   }

   QPSubproblem::~QPSubproblem() { }

   void QPSubproblem::initialize_statistics(Statistics& statistics, const Options& options) {
      if (this->use_regularization) {
         statistics.add_column("regulariz", Statistics::double_width - 4, options.get_int("statistics_regularization_column_order"));
      }
   }

   void QPSubproblem::generate_initial_iterate(const OptimizationProblem& problem, Iterate& initial_iterate) {
      if (this->enforce_linear_constraints_at_initial_iterate) {
         Preprocessing::enforce_linear_constraints(problem.model, initial_iterate.primals, initial_iterate.multipliers, *this->solver);
      }
   }

   void QPSubproblem::solve(Statistics& /*statistics*/, const OptimizationProblem& problem, Iterate& current_iterate,
         const Multipliers& current_multipliers, Direction& direction, const WarmstartInformation& warmstart_information) {
      const LagrangeNewtonSubproblem subproblem(problem, current_iterate, current_multipliers.constraints, *this->hessian_model, this->trust_region_radius);
      this->solver->solve_QP(subproblem, this->initial_point, direction, warmstart_information);
      InequalityConstrainedMethod::compute_dual_displacements(current_multipliers, direction.multipliers);
      this->number_subproblems_solved++;
      // reset the initial point
      this->initial_point.fill(0.);
   }
} // namespace