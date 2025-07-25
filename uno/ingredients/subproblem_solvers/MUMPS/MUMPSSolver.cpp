// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "MUMPSSolver.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "optimization/Direction.hpp"
#include "optimization/WarmstartInformation.hpp"
#if defined(HAS_MPI) && defined(MUMPS_PARALLEL)
#include "mpi.h"
#endif

#define USE_COMM_WORLD (-987654)

namespace uno {
   MUMPSSolver::MUMPSSolver(): DirectSymmetricIndefiniteLinearSolver() {
      this->workspace.sym = MUMPSSolver::GENERAL_SYMMETRIC;
#if defined(HAS_MPI) && defined(MUMPS_PARALLEL)
      // TODO load number of processes from option file
      this->mumps_structure.par = 1;
#else
      this->workspace.par = 1;
#endif
      this->workspace.job = MUMPSSolver::JOB_INIT;
      this->workspace.comm_fortran = USE_COMM_WORLD;
      dmumps_c(&this->workspace);
      // control parameters
      this->workspace.icntl[0] = -1;
      this->workspace.icntl[1] = -1;
      this->workspace.icntl[2] = -1;
      this->workspace.icntl[3] = 0;
      this->workspace.icntl[5] = 0; // no scaling
      this->workspace.icntl[7] = 0; // no scaling

      this->workspace.icntl[12] = 1;
      this->workspace.icntl[23] = 1; // ICNTL(24) controls the detection of “null pivot rows”

      /*
      // debug for MUMPS team
      this->mumps_structure.icntl[1] = 6; // ICNTL(2)=6
      this->mumps_structure.icntl[2] = 6; // ICNTL(3)=6
      this->mumps_structure.icntl[3] = 6; // ICNTL(4)=2
       */
   }

   MUMPSSolver::~MUMPSSolver() {
      this->workspace.job = MUMPSSolver::JOB_END;
      dmumps_c(&this->workspace);
   }

   void MUMPSSolver::initialize(const Subproblem& subproblem) {
      const size_t dimension = subproblem.number_variables + subproblem.number_constraints;

      // evaluations
      this->objective_gradient.resize(subproblem.number_variables);
      this->constraints.resize(subproblem.number_constraints);
      this->constraint_jacobian.resize(subproblem.number_constraints, subproblem.number_variables);

      // augmented system
      const size_t number_augmented_system_nonzeros = subproblem.number_augmented_system_nonzeros();
      const size_t regularization_size = subproblem.regularization_size();
      const size_t number_nonzeros = number_augmented_system_nonzeros + regularization_size;
      // reserve the COO sparse representation
      this->row_indices.reserve(number_nonzeros);
      this->column_indices.reserve(number_nonzeros);
      Vector<size_t> tmp_row_indices(number_nonzeros);
      Vector<size_t> tmp_column_indices(number_nonzeros);
      subproblem.compute_regularized_augmented_matrix_structure(tmp_row_indices, tmp_column_indices);
      // build vectors of int
      for (size_t nonzero_index: Range(number_nonzeros)) {
         this->row_indices[nonzero_index] = static_cast<int>(tmp_row_indices[nonzero_index]);
         this->column_indices[nonzero_index] = static_cast<int>(tmp_column_indices[nonzero_index]);
      }
      this->rhs.resize(dimension);
      this->solution.resize(dimension);

      // workspace
      this->workspace.n = static_cast<int>(dimension);
      this->workspace.nnz = static_cast<int>(number_nonzeros);
   }

   void MUMPSSolver::do_symbolic_analysis() {
      this->workspace.job = MUMPSSolver::JOB_ANALYSIS;
      // connect the local sparsity with the pointers in the structure
      this->workspace.irn = this->row_indices.data();
      this->workspace.jcn = this->column_indices.data();
      dmumps_c(&this->workspace);
      this->workspace.icntl[7] = 8; // ICNTL(8) = 8: recompute scaling before factorization
   }

   void MUMPSSolver::do_numerical_factorization(const Vector<double>& matrix_values) {
      this->workspace.job = MUMPSSolver::JOB_FACTORIZATION;
      this->workspace.a = const_cast<double*>(matrix_values.data());
      dmumps_c(&this->workspace);
   }

   void MUMPSSolver::solve_indefinite_system(const Vector<double>& /*matrix_values*/, const Vector<double>& rhs, Vector<double>& result) {
      result = rhs;
      this->workspace.rhs = result.data();
      this->workspace.job = MUMPSSolver::JOB_SOLVE;
      dmumps_c(&this->workspace);
   }

   void MUMPSSolver::solve_indefinite_system(Statistics& statistics, const Subproblem& subproblem, Direction& direction,
         const WarmstartInformation& warmstart_information) {
      // evaluate the functions at the current iterate
      if (warmstart_information.objective_changed) {
         subproblem.evaluate_objective_gradient(this->objective_gradient);
      }
      if (warmstart_information.constraints_changed) {
         subproblem.evaluate_constraints(this->constraints);
         // TODO subproblem.evaluate_jacobian(this->constraint_jacobian);
      }

      if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {
         // assemble the augmented matrix
         subproblem.assemble_augmented_matrix(statistics, this->matrix_values);
         // regularize the augmented matrix (this calls the analysis and the factorization)
         subproblem.regularize_augmented_matrix(statistics, this->matrix_values, subproblem.dual_regularization_factor(), *this);

         // assemble the RHS
         subproblem.assemble_augmented_rhs(this->objective_gradient, this->constraints, this->constraint_jacobian, this->rhs);
      }
      this->solve_indefinite_system(this->matrix_values, this->rhs, this->solution);
      // assemble the full primal-dual direction
      subproblem.assemble_primal_dual_direction(this->solution, direction);
      if (this->matrix_is_singular()) {
         direction.status = SubproblemStatus::INFEASIBLE;
      }
   }

   Inertia MUMPSSolver::get_inertia() const {
      const size_t number_negative_eigenvalues = this->number_negative_eigenvalues();
      const size_t number_zero_eigenvalues = this->number_zero_eigenvalues();
      const size_t number_positive_eigenvalues = static_cast<size_t>(this->workspace.n) - (number_negative_eigenvalues + number_zero_eigenvalues);
      return {number_positive_eigenvalues, number_negative_eigenvalues, number_zero_eigenvalues};
   }

   size_t MUMPSSolver::number_negative_eigenvalues() const {
      // INFOG(12)
      return static_cast<size_t>(this->workspace.infog[11]);
   }

   size_t MUMPSSolver::number_zero_eigenvalues() const {
      // INFOG(28)
      return static_cast<size_t>(this->workspace.infog[27]);
   }

   bool MUMPSSolver::matrix_is_singular() const {
      return (this->number_zero_eigenvalues() > 0);
   }

   size_t MUMPSSolver::rank() const {
      return this->workspace.n - this->number_zero_eigenvalues();
   }
} // namespace