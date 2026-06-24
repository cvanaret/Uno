// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HIGHSQUADRATICPROGRAM_H
#define UNO_HIGHSQUADRATICPROGRAM_H

#include <Highs.h>
#include "ingredients/subproblem_solvers/QuadraticProgram.hpp"
#include "ingredients/subproblem_solvers/SolverWorkspace.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declarations
   class Evaluations;
   class OptimizationProblem;
   class Statistics;
   class Subproblem;
   class WarmstartInformation;

   // HiGHS-native QuadraticProgram: the data is held directly in a HighsModel (CSC Jacobian + Hessian,
   // dense bounds and objective) inside the workspace. build() fills the model from the Subproblem.
   class HiGHSQuadraticProgram : public QuadraticProgram, public SolverWorkspace {
   public:
      HiGHSQuadraticProgram() = default;

      void initialize_memory(const Subproblem& subproblem) override;
      void fill(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) override;
      // data-driven build: dense objective gradient + COO constraint Jacobian + COO Lagrangian Hessian
      // (one triangle; empty for an LP). Converts COO to HiGHS' CSC layout internally.
      void fill(const Vector<double>& linear_objective,
         const Vector<uno_int>& jacobian_row_indices, const Vector<uno_int>& jacobian_column_indices,
         const Vector<double>& jacobian_values,
         const Vector<uno_int>& hessian_row_indices, const Vector<uno_int>& hessian_column_indices,
         const Vector<double>& hessian_values,
         const std::vector<double>& variables_lower_bounds, const std::vector<double>& variables_upper_bounds,
         const std::vector<double>& constraints_lower_bounds, const std::vector<double>& constraints_upper_bounds) override;
      [[nodiscard]] double compute_hessian_quadratic_form(const Subproblem& subproblem, const Vector<double>& vector) const override;

      HighsModel model;
      Vector<double> constraints{};
      Vector<double> linear_objective{};
      // constraint Jacobian in COO format
      Vector<uno_int> jacobian_row_indices{};
      Vector<uno_int> jacobian_column_indices{};
      Vector<double> jacobian_values{};
      Vector<size_t> jacobian_permutation_vector{};
      // Lagrangian Hessian in COO format
      Vector<uno_int> hessian_row_indices{};
      Vector<uno_int> hessian_column_indices{};
      Vector<double> hessian_values{};
      Vector<size_t> hessian_permutation_vector{};

   protected:
      void compute_jacobian_sparsity(const Subproblem& subproblem);
      void compute_hessian_sparsity(const Subproblem& subproblem);
      // data-driven setup: dense objective gradient + COO constraint Jacobian (row = constraint,
      // column = variable) + COO Lagrangian Hessian (one triangle; empty for an LP). Allocates the
      // HighsModel and converts the Jacobian/Hessian to HiGHS' CSC layout.
      void set_from_coo(size_t number_variables, size_t number_constraints, const Vector<double>& linear_objective,
         const Vector<uno_int>& jacobian_row_indices, const Vector<uno_int>& jacobian_column_indices,
         const Vector<double>& jacobian_values, const Vector<uno_int>& hessian_row_indices,
         const Vector<uno_int>& hessian_column_indices, const Vector<double>& hessian_values);

      void evaluate_functions(Statistics& statistics, const Subproblem& subproblem, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information);
      void evaluate_jacobian(const OptimizationProblem& problem, const Vector<double>& primals, Evaluations& evaluations);

      // build HiGHS' CSC Jacobian/Hessian from the COO arrays already stored in the *_row/column_indices members
      void build_csc_jacobian_from_coo(size_t number_variables, size_t number_constraints);
      void build_csc_hessian_from_coo(size_t number_variables);
      // scatter the COO values into the CSC value arrays using the sorting permutations
      void scatter_jacobian_values();
      void scatter_hessian_values();
   };
} // namespace

#endif // UNO_HIGHSQUADRATICPROGRAM_H
