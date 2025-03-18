// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "MUMPSSolver.hpp"
#include "ingredients/subproblems/LagrangeNewtonSubproblem.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "optimization/WarmstartInformation.hpp"

#if defined(HAS_MPI) && defined(MUMPS_PARALLEL)
#include "mpi.h"
#endif

#define USE_COMM_WORLD (-987654)

namespace uno {
   MUMPSSolver::MUMPSSolver(size_t number_variables, size_t number_constraints, size_t number_jacobian_nonzeros, size_t number_hessian_nonzeros) :
      DirectEqualityQPSolver<size_t, double>(),
      objective_gradient(number_variables),
      constraints(number_constraints),
      constraint_jacobian(number_constraints, number_variables), // TODO construct better
         hessian(number_variables, number_hessian_nonzeros, false, "COO"),
      dimension(number_variables + number_constraints),
      number_nonzeros(number_hessian_nonzeros + number_jacobian_nonzeros),
      augmented_matrix(this->dimension, this->number_nonzeros, true, "COO"),
      rhs(this->dimension),
      solution(this->dimension) {
      this->row_indices.reserve(number_nonzeros);
      this->column_indices.reserve(number_nonzeros);

      this->mumps_structure.sym = MUMPSSolver::GENERAL_SYMMETRIC;
#if defined(HAS_MPI) && defined(MUMPS_PARALLEL)
      // TODO load number of processes from option file
      this->mumps_structure.par = 1;
#else
      this->mumps_structure.par = 1;
#endif
      this->mumps_structure.job = MUMPSSolver::JOB_INIT;
      this->mumps_structure.comm_fortran = USE_COMM_WORLD;
      dmumps_c(&this->mumps_structure);
      // control parameters
      this->mumps_structure.icntl[0] = -1;
      this->mumps_structure.icntl[1] = -1;
      this->mumps_structure.icntl[2] = -1;
      this->mumps_structure.icntl[3] = 0;
      //this->mumps_structure.icntl[5] = 0; // no scaling
      //this->mumps_structure.icntl[7] = 0; // no scaling

      this->mumps_structure.icntl[12] = 1;
      this->mumps_structure.icntl[23] = 1; // ICNTL(24) controls the detection of “null pivot rows”

      /*
      // debug for MUMPS team
      this->mumps_structure.icntl[1] = 6; // ICNTL(2)=6
      this->mumps_structure.icntl[2] = 6; // ICNTL(3)=6
      this->mumps_structure.icntl[3] = 6; // ICNTL(4)=2
       */

      this->row_indices.reserve(number_nonzeros);
      this->column_indices.reserve(number_nonzeros);
   }

   MUMPSSolver::~MUMPSSolver() {
      this->mumps_structure.job = MUMPSSolver::JOB_END;
      dmumps_c(&this->mumps_structure);
   }

   void MUMPSSolver::do_symbolic_analysis(const SymmetricMatrix<size_t, double>& matrix) {
      this->mumps_structure.job = MUMPSSolver::JOB_ANALYSIS;
      this->mumps_structure.n = static_cast<int>(matrix.dimension());
      this->mumps_structure.nnz = static_cast<int>(matrix.number_nonzeros());
      this->mumps_structure.a = nullptr;
      this->save_sparsity_to_local_format(matrix);
      // connect the local sparsity with the pointers in the structure
      this->mumps_structure.irn = this->row_indices.data();
      this->mumps_structure.jcn = this->column_indices.data();
      this->mumps_structure.a = nullptr;
      dmumps_c(&this->mumps_structure);
      this->mumps_structure.icntl[7] = 8; // ICNTL(8) = 8: recompute scaling before factorization
   }

   void MUMPSSolver::do_numerical_factorization(const SymmetricMatrix<size_t, double>& matrix) {
      this->mumps_structure.job = MUMPSSolver::JOB_FACTORIZATION;
      this->mumps_structure.a = const_cast<double*>(matrix.data_pointer());
      dmumps_c(&this->mumps_structure);
   }

   SubproblemStatus MUMPSSolver::solve_equality_constrained_QP(Statistics& statistics, LagrangeNewtonSubproblem& subproblem,
         const Vector<double>& /*initial_point*/, Vector<double>& direction_primals, Multipliers& direction_multipliers, double& /*subproblem_objective*/,
         WarmstartInformation& warmstart_information) {
      // set up the augmented system
      subproblem.assemble_augmented_matrix(statistics, this->objective_gradient, this->constraints, this->constraint_jacobian, this->hessian,
         this->augmented_matrix, *this, warmstart_information);
      subproblem.assemble_augmented_rhs(this->objective_gradient, this->constraints, this->constraint_jacobian, this->rhs, warmstart_information);
      // solve the augmented system
      this->solve_indefinite_linear_system();
      // form the primal-dual direction (note the minus sign for the multipliers)
      direction_primals = view(this->solution, 0, subproblem.number_variables);
      direction_multipliers.constraints = view(-this->solution, subproblem.number_variables, subproblem.number_variables + subproblem.number_constraints);
      return SubproblemStatus::OPTIMAL; // TODO
   }

   void MUMPSSolver::solve_indefinite_linear_system() {
      // copy rhs into solution (overwritten by MUMPS)
      this->solution = this->rhs;
      this->mumps_structure.rhs = this->solution.data();
      this->mumps_structure.job = MUMPSSolver::JOB_SOLVE;
      dmumps_c(&this->mumps_structure);
   }

   std::tuple<size_t, size_t, size_t> MUMPSSolver::get_inertia() const {
      const size_t number_negative_eigenvalues = this->number_negative_eigenvalues();
      const size_t number_zero_eigenvalues = this->number_zero_eigenvalues();
      const size_t number_positive_eigenvalues =
            static_cast<size_t>(this->mumps_structure.n) - (number_negative_eigenvalues + number_zero_eigenvalues);
      return std::make_tuple(number_positive_eigenvalues, number_negative_eigenvalues, number_zero_eigenvalues);
   }

   size_t MUMPSSolver::number_negative_eigenvalues() const {
      // INFOG(12)
      return static_cast<size_t>(this->mumps_structure.infog[11]);
   }

   size_t MUMPSSolver::number_zero_eigenvalues() const {
      // INFOG(28)
      return static_cast<size_t>(this->mumps_structure.infog[27]);
   }

   bool MUMPSSolver::matrix_is_singular() const {
      return (this->number_zero_eigenvalues() > 0);
   }

   size_t MUMPSSolver::rank() const {
      return this->dimension - this->number_zero_eigenvalues();
   }

   void MUMPSSolver::save_sparsity_to_local_format(const SymmetricMatrix<size_t, double>& matrix) {
      // build the internal matrix representation
      this->row_indices.clear();
      this->column_indices.clear();
      for (const auto[row_index, column_index, _]: matrix) {
         this->row_indices.emplace_back(static_cast<int>(row_index + this->fortran_shift));
         this->column_indices.emplace_back(static_cast<int>(column_index + this->fortran_shift));
      }
   }
} // namespace
