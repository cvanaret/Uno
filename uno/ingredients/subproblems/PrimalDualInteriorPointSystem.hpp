// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALDUALINTERIORPOINTSYSTEM_H
#define UNO_PRIMALDUALINTERIORPOINTSYSTEM_H

#include "LagrangeNewtonSubproblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "reformulation/OptimizationProblem.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/Options.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   // forward declaration
   template <typename ElementType>
   class Vector;

   class PrimalDualInteriorPointSystem: public LagrangeNewtonSubproblem {
   public:
      PrimalDualInteriorPointSystem(const OptimizationProblem& problem, Iterate& current_iterate, const Multipliers& current_multipliers,
            double barrier_parameter, bool use_regularization, double trust_region_radius, const Options& options);

      void evaluate_matrix(SymmetricMatrix<size_t, double>& matrix, const WarmstartInformation& warmstart_information) const;
      void evaluate_right_hand_side(Vector<double>& rhs, const WarmstartInformation& warmstart_information) const;

      const size_t number_variables;

   protected:
      const double barrier_parameter;
      const double damping_factor{1e-4}; // TODO
   };
} // namespace

#endif // UNO_PRIMALDUALINTERIORPOINTSYSTEM_H