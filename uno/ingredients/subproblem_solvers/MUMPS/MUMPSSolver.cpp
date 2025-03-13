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
         DirectSymmetricIndefiniteLinearSolver<size_t, double>(number_variables + number_constraints),
         objective_gradient(number_variables),
         constraints(number_constraints),
         constraint_jacobian(number_constraints, number_variables), // TODO construct better
         hessian(number_variables, number_hessian_nonzeros, false, "COO"),
         dimension(number_variables + number_constraints),
         number_nonzeros(number_hessian_nonzeros + number_jacobian_nonzeros),
         augmented_matrix(this->dimension, this->number_nonzeros, true, "COO"),
         rhs(this->dimension) {
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

   void MUMPSSolver::solve_indefinite_system(const SymmetricMatrix<size_t, double>& /*matrix*/, const Vector<double>& rhs, Vector<double>& result) {
      result = rhs;
      this->mumps_structure.rhs = result.data();
      this->mumps_structure.job = MUMPSSolver::JOB_SOLVE;
      dmumps_c(&this->mumps_structure);
   }

   void MUMPSSolver::solve_indefinite_system(Statistics& statistics, LagrangeNewtonSubproblem& subproblem, Vector<double>& result,
         WarmstartInformation& warmstart_information) {
      this->set_up_subproblem(statistics, subproblem, warmstart_information);
      this->solve_indefinite_system(this->augmented_matrix, this->rhs, result);
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

   void MUMPSSolver::set_up_subproblem(Statistics& statistics, LagrangeNewtonSubproblem& subproblem, WarmstartInformation& warmstart_information) {
      // objective gradient
      if (warmstart_information.objective_changed) {
         subproblem.evaluate_objective_gradient(this->objective_gradient);
      }

      // constraints and Jacobian
      if (warmstart_information.constraints_changed) {
         subproblem.evaluate_constraints(this->constraints);
      }
      if (warmstart_information.constraint_jacobian_changed) {
         subproblem.evaluate_constraint_jacobian(this->constraint_jacobian);
      }

      // Lagrangian Hessian
      if (warmstart_information.objective_changed || warmstart_information.constraints_changed || warmstart_information.constraint_jacobian_changed) {
         DEBUG << "Evaluating problem Hessian\n";
         this->hessian.reset();
         subproblem.evaluate_hessian(this->hessian);
      }

      if (warmstart_information.objective_changed || warmstart_information.constraint_jacobian_changed) {
         // form the KKT matrix
         this->augmented_matrix.set_dimension(subproblem.number_variables + subproblem.number_constraints);
         this->augmented_matrix.reset();
         // copy the Lagrangian Hessian in the top left block
         for (const auto [row_index, column_index, element]: this->hessian) {
            this->augmented_matrix.insert(element, row_index, column_index);
         }

         // Jacobian of general constraints
         for (size_t column_index: Range(subproblem.number_constraints)) {
            for (const auto [row_index, derivative]: this->constraint_jacobian[column_index]) {
               this->augmented_matrix.insert(derivative, row_index, subproblem.number_variables + column_index);
            }
            this->augmented_matrix.finalize_column(column_index);
         }
      }

      // possibly assemble augmented system and perform analysis
      if (warmstart_information.hessian_sparsity_changed || warmstart_information.jacobian_sparsity_changed) {
         DEBUG << "Augmented matrix:\n" << this->augmented_matrix;
         DEBUG << "Performing symbolic analysis of the augmented matrix\n";
         this->do_symbolic_analysis(this->augmented_matrix);
         warmstart_information.hessian_sparsity_changed = warmstart_information.jacobian_sparsity_changed = false;
      }
      if (warmstart_information.objective_changed || warmstart_information.constraint_jacobian_changed) {
         DEBUG << "Performing numerical factorization of the augmented matrix\n";
         this->do_numerical_factorization(this->augmented_matrix);
         subproblem.regularize_matrix(statistics, *this, this->augmented_matrix);
      }
      this->assemble_augmented_rhs(subproblem); // TODO add conditions
   }

   void MUMPSSolver::assemble_augmented_rhs(LagrangeNewtonSubproblem& subproblem) {
      // Lagrangian gradient
      subproblem.compute_lagrangian_gradient(this->objective_gradient, this->constraint_jacobian, this->rhs);
      for (size_t variable_index: Range(subproblem.number_variables)) {
         this->rhs[variable_index] = -this->rhs[variable_index];
      }

      // constraints
      for (size_t constraint_index: Range(subproblem.number_constraints)) {
         this->rhs[subproblem.number_variables + constraint_index] = -this->constraints[constraint_index];
      }
      DEBUG2 << "RHS: "; print_vector(DEBUG2, view(this->rhs, 0, subproblem.number_variables + subproblem.number_constraints));
      DEBUG << '\n';
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
