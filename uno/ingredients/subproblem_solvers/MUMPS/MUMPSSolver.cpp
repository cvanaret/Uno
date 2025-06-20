// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "MUMPSSolver.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/WarmstartInformation.hpp"
#if defined(HAS_MPI) && defined(MUMPS_PARALLEL)
#include "mpi.h"
#endif

#define USE_COMM_WORLD (-987654)

namespace uno {
   MUMPSSolver::MUMPSSolver(): DirectSymmetricIndefiniteLinearSolver() {
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
   }

   MUMPSSolver::~MUMPSSolver() {
      this->mumps_structure.job = MUMPSSolver::JOB_END;
      dmumps_c(&this->mumps_structure);
   }

   void MUMPSSolver::initialize_memory(size_t number_variables, size_t number_constraints, size_t number_hessian_nonzeros,
         size_t regularization_size) {
      const size_t dimension = number_variables + number_constraints;
      this->dimension = dimension;

      // reserve the COO sparse representation
      const size_t number_nonzeros = number_hessian_nonzeros + regularization_size;
      this->row_indices.reserve(number_nonzeros);
      this->column_indices.reserve(number_nonzeros);

      // evaluations
      this->objective_gradient.reserve(number_variables);
      this->constraints.resize(number_constraints);
      this->constraint_jacobian.resize(number_constraints, number_variables);
      this->augmented_matrix = SymmetricMatrix<size_t, double>("COO", dimension, number_hessian_nonzeros, regularization_size);
      this->rhs.resize(dimension);
   }

   void MUMPSSolver::do_symbolic_analysis(const SymmetricMatrix<size_t, double>& matrix) {
      this->mumps_structure.job = MUMPSSolver::JOB_ANALYSIS;
      this->mumps_structure.n = static_cast<int>(matrix.dimension());
      this->mumps_structure.nnz = static_cast<int>(matrix.number_nonzeros());
      this->save_sparsity_to_local_format(matrix);
      // connect the local sparsity with the pointers in the structure
      this->mumps_structure.irn = this->row_indices.data();
      this->mumps_structure.jcn = this->column_indices.data();
      //this->mumps_structure.a = const_cast<double*>(matrix.data_pointer());
      dmumps_c(&this->mumps_structure);
      this->mumps_structure.icntl[7] = 8; // ICNTL(8) = 8: recompute scaling before factorization
   }

   void MUMPSSolver::do_numerical_factorization(const SymmetricMatrix<size_t, double>& matrix) {
      this->mumps_structure.job = MUMPSSolver::JOB_FACTORIZATION;
      this->mumps_structure.a = const_cast<double*>(matrix.data_pointer());
      dmumps_c(&this->mumps_structure);
   }

   void MUMPSSolver::solve_indefinite_system(const SymmetricMatrix<size_t, double>& /*matrix*/, const Vector<double>& rhs, Vector<double>& result) {
      result = rhs;
      this->mumps_structure.rhs = result.data();
      this->mumps_structure.job = MUMPSSolver::JOB_SOLVE;
      dmumps_c(&this->mumps_structure);
   }

   void MUMPSSolver::solve_indefinite_system(Statistics& statistics, const Subproblem& subproblem, Vector<double>& solution,
         const WarmstartInformation& warmstart_information) {
      // evaluate the functions at the current iterate
      if (warmstart_information.objective_changed) {
         subproblem.evaluate_objective_gradient(this->objective_gradient);
      }
      if (warmstart_information.constraints_changed) {
         subproblem.evaluate_constraints(this->constraints);
         subproblem.evaluate_jacobian(this->constraint_jacobian);
      }

      if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {
         // assemble the augmented matrix
         this->augmented_matrix.reset();
         subproblem.assemble_augmented_matrix(statistics, this->augmented_matrix, this->constraint_jacobian);
         // regularize the augmented matrix (this calls the analysis and the factorization)
         subproblem.regularize_augmented_matrix(statistics, this->augmented_matrix, subproblem.dual_regularization_factor(), *this);

         // assemble the RHS
         subproblem.assemble_augmented_rhs(this->objective_gradient, this->constraints, this->constraint_jacobian, this->rhs);
      }
      this->solve_indefinite_system(this->augmented_matrix, this->rhs, solution);
   }

   Inertia MUMPSSolver::get_inertia() const {
      const size_t number_negative_eigenvalues = this->number_negative_eigenvalues();
      const size_t number_zero_eigenvalues = this->number_zero_eigenvalues();
      const size_t number_positive_eigenvalues = static_cast<size_t>(this->mumps_structure.n) - (number_negative_eigenvalues + number_zero_eigenvalues);
      return {number_positive_eigenvalues, number_negative_eigenvalues, number_zero_eigenvalues};
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
      for (const auto [row_index, column_index, _]: matrix) {
         this->row_indices.emplace_back(static_cast<int>(row_index + this->fortran_shift));
         this->column_indices.emplace_back(static_cast<int>(column_index + this->fortran_shift));
      }
   }
} // namespace