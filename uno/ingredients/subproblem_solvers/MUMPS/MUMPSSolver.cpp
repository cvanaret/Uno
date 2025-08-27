// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "MUMPSSolver.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#if defined(HAS_MPI) && defined(MUMPS_PARALLEL)
#include "mpi.h"
#endif

#define USE_COMM_WORLD (-987654)

namespace uno {
   MUMPSSolver::MUMPSSolver(): DirectSymmetricIndefiniteLinearSolver() {
      this->workspace.sym = MUMPSSolver::GENERAL_SYMMETRIC;
#if defined(HAS_MPI) && defined(MUMPS_PARALLEL)
      // TODO load number of processes from option file
      this->workspace.par = 1;
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
      this->workspace.icntl[1] = 6; // ICNTL(2)=6
      this->workspace.icntl[2] = 6; // ICNTL(3)=6
      this->workspace.icntl[3] = 6; // ICNTL(4)=2
       */
   }

   MUMPSSolver::~MUMPSSolver() {
      this->workspace.job = MUMPSSolver::JOB_END;
      dmumps_c(&this->workspace);
   }

   void MUMPSSolver::initialize_hessian(const Subproblem& subproblem) {
      this->evaluation_space.initialize_hessian(subproblem);

      // workspace*
      const size_t dimension = subproblem.number_variables;
      this->workspace.n = static_cast<int>(dimension);
      this->workspace.nnz = static_cast<int>(this->evaluation_space.number_matrix_nonzeros);
   }

   void MUMPSSolver::initialize_augmented_system(const Subproblem& subproblem) {
      this->evaluation_space.initialize_augmented_system(subproblem);

      // workspace
      const size_t dimension = subproblem.number_variables + subproblem.number_constraints;
      this->workspace.n = static_cast<int>(dimension);
      this->workspace.nnz = static_cast<int>(this->evaluation_space.number_matrix_nonzeros);
   }

   void MUMPSSolver::do_symbolic_analysis() {
      assert(!this->analysis_performed);

      this->workspace.job = MUMPSSolver::JOB_ANALYSIS;
      // connect the local sparsity with the pointers in the workspace
      this->workspace.irn = this->evaluation_space.matrix_row_indices.data();
      this->workspace.jcn = this->evaluation_space.matrix_column_indices.data();
      dmumps_c(&this->workspace);
      this->workspace.icntl[7] = 8; // ICNTL(8) = 8: recompute scaling before factorization
      this->analysis_performed = true;
   }

   void MUMPSSolver::do_numerical_factorization(const double* matrix_values) {
      assert(this->analysis_performed);

      this->workspace.job = MUMPSSolver::JOB_FACTORIZATION;
      this->workspace.a = const_cast<double*>(matrix_values);
      dmumps_c(&this->workspace);
      this->factorization_performed = true;
   }

   void MUMPSSolver::solve_indefinite_system(const Vector<double>& /*matrix_values*/, const Vector<double>& rhs, Vector<double>& result) {
      assert(this->factorization_performed);

      result = rhs;
      this->workspace.rhs = result.data();
      this->workspace.job = MUMPSSolver::JOB_SOLVE;
      dmumps_c(&this->workspace);
   }

   void MUMPSSolver::solve_indefinite_system(Statistics& statistics, const Subproblem& subproblem, Direction& direction,
         const WarmstartInformation& warmstart_information) {
      // set up the linear system by evaluating the functions at the current iterate
      this->evaluation_space.set_up_linear_system(statistics, subproblem, *this, warmstart_information);
      // solve the linear system
      this->solve_indefinite_system(this->evaluation_space.matrix_values, this->evaluation_space.rhs, this->evaluation_space.solution);
      // assemble the full primal-dual direction
      subproblem.assemble_primal_dual_direction(this->evaluation_space.solution, direction);
      if (this->matrix_is_singular()) {
         direction.status = SubproblemStatus::INFEASIBLE;
      }
   }

   Inertia MUMPSSolver::get_inertia() const {
      // rank = number_positive_eigenvalues + number_negative_eigenvalues
      // n = rank + number_zero_eigenvalues
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

   EvaluationSpace& MUMPSSolver::get_evaluation_space() {
      return this->evaluation_space;
   }
} // namespace