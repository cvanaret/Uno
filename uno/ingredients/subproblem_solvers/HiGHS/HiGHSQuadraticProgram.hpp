// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HIGHSQUADRATICPROGRAM_H
#define UNO_HIGHSQUADRATICPROGRAM_H

#include <vector>
#include <interfaces/highs_c_api.h>
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

   // HiGHS-native QuadraticProgram: the data is held directly as the raw arrays that Highs_passModel
   // consumes (dense objective and bounds, CSC constraint Jacobian, triangular CSC Hessian), all owned
   // here inside the workspace. build()/fill() populate them from the Subproblem or from raw data.
   class HiGHSQuadraticProgram : public QuadraticProgram, public SolverWorkspace {
   public:
      HiGHSQuadraticProgram() = default;

      void initialize_memory(const Subproblem& subproblem) override;
      void fill(Statistics& statistics, const Subproblem& subproblem, const Iterate& current_iterate, double trust_region_radius,
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
      [[nodiscard]] double compute_hessian_quadratic_form(const Subproblem& subproblem, const Iterate& current_iterate,
         const Vector<double>& vector) const override;

      // HiGHS C-API model data (owned here; handed to Highs_passModel each solve).
      // Dimensions are number_variables (columns) and number_constraints (rows).
      std::vector<double> col_cost{};    // dense objective gradient
      std::vector<double> col_lower{};   // variable lower bounds
      std::vector<double> col_upper{};   // variable upper bounds
      std::vector<double> row_lower{};   // constraint lower bounds
      std::vector<double> row_upper{};   // constraint upper bounds
      // column-wise (CSC) constraint Jacobian
      std::vector<HighsInt> a_start{};   // column starts, length number_variables + 1
      std::vector<HighsInt> a_index{};   // row (constraint) indices
      std::vector<double> a_value{};
      // lower-triangular CSC Lagrangian Hessian; left empty for an LP
      std::vector<HighsInt> q_start{};   // column starts, length number_variables + 1
      std::vector<HighsInt> q_index{};   // row indices
      std::vector<double> q_value{};
      // solution buffers filled by Highs_getSolution (preallocated in initialize_memory / fill)
      std::vector<double> col_value{};
      std::vector<double> col_dual{};
      std::vector<double> row_value{};
      std::vector<double> row_dual{};

      Vector<double> constraints{};
      Vector<double> linear_objective{};
      // constraint Jacobian in COO format
      // Vector<uno_int> jacobian_row_indices{};
      // Vector<uno_int> jacobian_column_indices{};
      Vector<double> jacobian_values{};
      Vector<size_t> jacobian_permutation_vector{};
      // Lagrangian Hessian in COO format
      Vector<uno_int> hessian_row_indices{};
      Vector<uno_int> hessian_column_indices{};
      Vector<double> hessian_values{}; // workspace
      Vector<size_t> hessian_permutation_vector{};

   protected:
      void compute_jacobian_sparsity(const Subproblem& subproblem);
      void compute_hessian_sparsity(const Subproblem& subproblem);
      // data-driven setup: dense objective gradient + COO constraint Jacobian (row = constraint,
      // column = variable) + COO Lagrangian Hessian (one triangle; empty for an LP). Sizes the arrays
      // and converts the Jacobian/Hessian to HiGHS' CSC layout.
      void set_from_coo(size_t number_variables, size_t number_constraints, const Vector<double>& linear_objective,
         const Vector<uno_int>& jacobian_row_indices, const Vector<uno_int>& jacobian_column_indices,
         const Vector<double>& jacobian_values, const Vector<uno_int>& hessian_row_indices,
         const Vector<uno_int>& hessian_column_indices, const Vector<double>& hessian_values);

      void evaluate_functions(Statistics& statistics, const Subproblem& subproblem, const Iterate& current_iterate,
         Evaluations& current_evaluations, const WarmstartInformation& warmstart_information);
      void evaluate_jacobian(const OptimizationProblem& problem, const Vector<double>& primals, Evaluations& evaluations);

      // build HiGHS' CSC Jacobian/Hessian from the COO arrays already stored in the *_row/column_indices members
      void build_csc_jacobian_from_coo(size_t number_variables, size_t number_constraints,
         const Vector<uno_int>& jacobian_row_indices, const Vector<uno_int>& jacobian_column_indices);
      void build_csc_hessian_from_coo(size_t number_variables);
      // scatter the COO values into the CSC value arrays using the sorting permutations
      void scatter_jacobian_values();
      void scatter_hessian_values();
   };
} // namespace

#endif // UNO_HIGHSQUADRATICPROGRAM_H
