// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "MA57Solver.hpp"
#include "ingredients/subproblems/LagrangeNewtonSubproblem.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/Logger.hpp"
#include "fortran_interface.h"

#define MA57ID FC_GLOBAL(ma57id, MA57ID)
#define MA57AD FC_GLOBAL(ma57ad, MA57AD)
#define MA57BD FC_GLOBAL(ma57bd, MA57BD)
#define MA57CD FC_GLOBAL(ma57cd, MA57CD)
#define MA57DD FC_GLOBAL(ma57dd, MA57DD)

namespace uno {
   extern "C" {
      // default values of controlling parameters
      void MA57ID(double cntl[], int icntl[]);
      // symbolic analysis
      void MA57AD(const int* n, const int* ne, const int irn[], const int jcn[], const int* lkeep, int keep[], int iwork[], int icntl[], int info[],
            double rinfo[]);
      // numerical factorization
      void MA57BD(const int* n, int* ne, const double a[], /* out */ double fact[], const int* lfact, /* out */ int ifact[], const int* lifact,
            const int* lkeep, const int keep[], int iwork[], int icntl[], double cntl[], /* out */ int info[], /* out */ double rinfo[]);
      // linear system solve without iterative refinement
      void MA57CD(const int* job, const int* n, double fact[], int* lfact, int ifact[], int* lifact, const int* nrhs, double rhs[], const int* lrhs,
            double work[], int* lwork, int iwork[], int icntl[], int info[]);
      // linear system solve with iterative refinement
      void MA57DD(const int* job, const int* n, int* ne, const double a[], const int irn[], const int jcn[], double fact[], int* lfact, int ifact[],
            int* lifact, const double rhs[], double x[], double resid[], double work[], int iwork[], int icntl[], double cntl[], int info[],
            double rinfo[]);
   }

   MA57Solver::MA57Solver(size_t number_variables, size_t number_constraints, size_t number_jacobian_nonzeros, size_t number_hessian_nonzeros):
            DirectSymmetricIndefiniteLinearSolver<size_t, double>(number_variables + number_constraints),
            objective_gradient(number_variables),
            constraints(number_constraints),
            constraint_jacobian(number_constraints, number_variables), // TODO construct better
            hessian(number_variables, number_hessian_nonzeros, false, "COO"),
            dimension(number_variables + number_constraints),
            number_nonzeros(number_hessian_nonzeros + number_jacobian_nonzeros),
            augmented_matrix(this->dimension, this->number_nonzeros, true, "COO"),
            rhs(this->dimension),
            lkeep(static_cast<int>(5 * this->dimension + this->number_nonzeros + std::max(this->dimension, this->number_nonzeros) + 42)),
            keep(static_cast<size_t>(lkeep)),
            iwork(5 * this->dimension),
            lwork(static_cast<int>(1.2 * static_cast<double>(this->dimension))),
            work(static_cast<size_t>(this->lwork)), residuals(this->dimension) {
      this->row_indices.reserve(this->number_nonzeros);
      this->column_indices.reserve(this->number_nonzeros);
      // set the default values of the controlling parameters
      MA57ID(this->cntl.data(), this->icntl.data());
      // suppress warning messages
      this->icntl[4] = 0;
      // iterative refinement enabled
      this->icntl[8] = 1;
   }

   void MA57Solver::do_symbolic_analysis(const SymmetricMatrix<size_t, double>& matrix) {
      assert(matrix.dimension() <= this->dimension && "MA57Solver: the dimension of the matrix is larger than the preallocated size");
      assert(matrix.number_nonzeros() <= this->row_indices.capacity() &&
             "MA57Solver: the number of nonzeros of the matrix is larger than the preallocated size");

      // build the internal matrix representation
      this->save_sparsity_pattern_internally(matrix);

      const int n = static_cast<int>(matrix.dimension());
      const int nnz = static_cast<int>(matrix.number_nonzeros());

      // symbolic analysis
      MA57AD(/* const */ &n,
            /* const */ &nnz,
            /* const */ this->row_indices.data(),
            /* const */ this->column_indices.data(),
            /* const */ &this->lkeep,
            /* const */ this->keep.data(),
            /* out */ this->iwork.data(),
            /* const */ this->icntl.data(),
            /* out */ this->info.data(),
            /* out */ this->rinfo.data());

      assert(0 <= this->info[0] && "MA57: the symbolic analysis failed");
      if (0 < this->info[0]) {
         WARNING << "MA57 has issued a warning: info(1) = " << this->info[0] << '\n';
      }

      // get LFACT and LIFACT and resize FACT and IFACT (no effect if resized to <= size)
      int lfact = 2 * this->info[8];
      int lifact = 2 * this->info[9];
      this->fact.resize(static_cast<size_t>(lfact));
      this->ifact.resize(static_cast<size_t>(lifact));

      // store the sizes of the symbolic analysis
      this->factorization = {n, nnz, lfact, lifact};
   }

   void MA57Solver::do_numerical_factorization(const SymmetricMatrix<size_t, double>& matrix) {
      assert(matrix.dimension() <= this->dimension && "MA57Solver: the dimension of the matrix is larger than the preallocated size");
      assert(this->factorization.nnz == static_cast<int>(matrix.number_nonzeros()) && "MA57Solver: the numbers of nonzeros do not match");

      const int n = static_cast<int>(matrix.dimension());
      int nnz = static_cast<int>(matrix.number_nonzeros());

      // numerical factorization
      MA57BD(&n,
            &nnz,
            /* const */ matrix.data_pointer(),
            /* out */ this->fact.data(),
            /* const */ &this->factorization.lfact,
            /* out */ this->ifact.data(),
            /* const */ &this->factorization.lifact,
            /* const */ &this->lkeep,
            /* const */ this->keep.data(), this->iwork.data(), this->icntl.data(), this->cntl.data(),
            /* out */ this->info.data(),
            /* out */ this->rinfo.data());
   }

   void MA57Solver::solve_indefinite_system(const SymmetricMatrix<size_t, double>& matrix, const Vector<double>& rhs, Vector<double>& result) {
      const int n = static_cast<int>(matrix.dimension());
      int nnz = static_cast<int>(matrix.number_nonzeros());
      const int lrhs = n; // integer, length of rhs

      // solve the linear system
      if (this->use_iterative_refinement) {
         MA57DD(&this->job, &n, &nnz, matrix.data_pointer(), this->row_indices.data(), this->column_indices.data(),
               this->fact.data(), &this->factorization.lfact, this->ifact.data(), &this->factorization.lifact,
               rhs.data(), result.data(), this->residuals.data(), this->work.data(), this->iwork.data(), this->icntl.data(),
               this->cntl.data(), this->info.data(), this->rinfo.data());
      }
      else {
         // copy rhs into result (overwritten by MA57)
         result = rhs;

         MA57CD(&this->job, &n, this->fact.data(), &this->factorization.lfact, this->ifact.data(),
               &this->factorization.lifact, &this->nrhs, result.data(), &lrhs, this->work.data(), &this->lwork, this->iwork.data(),
               this->icntl.data(), this->info.data());
      }
   }

    void MA57Solver::solve_indefinite_system(Statistics& statistics, LagrangeNewtonSubproblem& subproblem, Vector<double>& result,
         WarmstartInformation& warmstart_information) {
       this->set_up_subproblem(statistics, subproblem, warmstart_information);
       this->solve_indefinite_system(this->augmented_matrix, this->rhs, result);
    }

   void MA57Solver::set_up_subproblem(Statistics& statistics, LagrangeNewtonSubproblem& subproblem, WarmstartInformation& warmstart_information) {
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

   void MA57Solver::assemble_augmented_rhs(LagrangeNewtonSubproblem& subproblem) {
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

   std::tuple<size_t, size_t, size_t> MA57Solver::get_inertia() const {
      // rank = number_positive_eigenvalues + number_negative_eigenvalues
      // n = rank + number_zero_eigenvalues
      const size_t rank = this->rank();
      const size_t number_negative_eigenvalues = this->number_negative_eigenvalues();
      const size_t number_positive_eigenvalues = rank - number_negative_eigenvalues;
      const size_t number_zero_eigenvalues = static_cast<size_t>(this->factorization.n) - rank;
      return std::make_tuple(number_positive_eigenvalues, number_negative_eigenvalues, number_zero_eigenvalues);
   }

   size_t MA57Solver::number_negative_eigenvalues() const {
      return static_cast<size_t>(this->info[23]);
   }

   /*
   bool MA57Solver::matrix_is_positive_definite() const {
      // positive definite = non-singular and no negative eigenvalues
      return not this->matrix_is_singular() && this->number_negative_eigenvalues() == 0;
   }
   */

   bool MA57Solver::matrix_is_singular() const {
      return (this->info[0] == 4);
   }

   size_t MA57Solver::rank() const {
      return static_cast<size_t>(this->info[24]);
   }

   void MA57Solver::save_sparsity_pattern_internally(const SymmetricMatrix<size_t, double>& matrix) {
      // build the internal matrix representation
      this->row_indices.clear();
      this->column_indices.clear();
      for (const auto [row_index, column_index, _]: matrix) {
         this->row_indices.emplace_back(static_cast<int>(row_index + this->fortran_shift));
         this->column_indices.emplace_back(static_cast<int>(column_index + this->fortran_shift));
      }
   }
} // namespace