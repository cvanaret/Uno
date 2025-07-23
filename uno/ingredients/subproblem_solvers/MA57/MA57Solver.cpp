// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "MA57Solver.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/Logger.hpp"
#include "fortran_interface.h"

#define MA57_set_default_parameters FC_GLOBAL(ma57id, MA57ID)
#define MA57_symbolic_analysis FC_GLOBAL(ma57ad, MA57AD)
#define MA57_numerical_factorization FC_GLOBAL(ma57bd, MA57BD)
#define MA57_linear_solve FC_GLOBAL(ma57cd, MA57CD)
#define MA57_linear_solve_with_iterative_refinement FC_GLOBAL(ma57dd, MA57DD)

namespace uno {
   extern "C" {
   // MA57
   // default values of controlling parameters
   void MA57_set_default_parameters(double cntl[], int icntl[]);

   // symbolic analysis
   void MA57_symbolic_analysis(const int* n, const int* ne, const int irn[], const int jcn[], const int* lkeep,
      /*const*/ int keep[], int iwork[], /* const */ int icntl[], int info[], double rinfo[]);

   // numerical factorization
   void MA57_numerical_factorization(const int* n, int* ne, const double a[], double fact[], const int* lfact,
      int ifact[], const int* lifact, const int* lkeep, const int keep[], int iwork[], int icntl[], double cntl[],
      int info[], double rinfo[]);

   // linear system solve without iterative refinement
   void MA57_linear_solve(const int* job, const int* n, double fact[], int* lfact, int ifact[], int* lifact, const int* nrhs,
      double rhs[], const int* lrhs, double work[], int* lwork, int iwork[], int icntl[], int info[]);

   // linear system solve with iterative refinement
   void MA57_linear_solve_with_iterative_refinement(const int* job, const int* n, int* ne, const double a[], const int irn[],
      const int jcn[], double fact[], int* lfact, int ifact[], int* lifact, const double rhs[], double x[], double resid[],
      double work[], int iwork[], int icntl[], double cntl[], int info[], double rinfo[]);
   }

   MA57Solver::MA57Solver(): DirectSymmetricIndefiniteLinearSolver() {
      // set the default values of the controlling parameters
      MA57_set_default_parameters(this->workspace.cntl.data(), this->workspace.icntl.data());
      // suppress warning messages
      this->workspace.icntl[4] = 0;
      // iterative refinement enabled
      this->workspace.icntl[8] = 1;
   }

   void MA57Solver::initialize_memory(size_t number_variables, size_t number_constraints, size_t number_hessian_nonzeros,
         size_t regularization_size) {
      const size_t dimension = number_variables + number_constraints;
      const size_t number_nonzeros = number_hessian_nonzeros + regularization_size;
      this->dimension = dimension;

      // reserve the COO sparse representation
      this->row_indices.reserve(number_nonzeros);
      this->column_indices.reserve(number_nonzeros);

      // evaluations
      this->objective_gradient.resize(number_variables);
      this->constraints.resize(number_constraints);
      this->constraint_jacobian.resize(number_constraints, number_variables);

      // augmented system
      this->augmented_matrix = SparseSymmetricMatrix<COOFormat<size_t, double>>(dimension, number_hessian_nonzeros, regularization_size);
      this->rhs.resize(dimension);
      this->solution.resize(dimension);

      // workspace
      this->workspace.n = static_cast<int>(dimension);
      this->workspace.nnz = static_cast<int>(number_nonzeros);
      this->workspace.lkeep = static_cast<int>(5 * dimension + number_nonzeros + std::max(dimension, number_nonzeros) + 42);
      this->workspace.keep.resize(static_cast<size_t>(this->workspace.lkeep));
      this->workspace.iwork.resize(5 * dimension);
      this->workspace.lwork = static_cast<int>(1.2 * static_cast<double>(dimension));
      this->workspace.work.resize(static_cast<size_t>(this->workspace.lwork));
      this->workspace.residuals.resize(dimension);
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
      MA57_symbolic_analysis(&n, &nnz, this->row_indices.data(), this->column_indices.data(), &this->workspace.lkeep,
         this->workspace.keep.data(), this->workspace.iwork.data(), this->workspace.icntl.data(), this->workspace.info.data(),
         this->workspace.rinfo.data());

      assert(0 <= this->workspace.info[0] && "MA57: the symbolic analysis failed");
      if (0 < this->workspace.info[0]) {
         WARNING << "MA57 has issued a warning: info(1) = " << workspace.info[0] << '\n';
      }

      // get LFACT and LIFACT and resize FACT and IFACT (no effect if resized to <= size)
      int lfact = 2 * this->workspace.info[8];
      int lifact = 2 * this->workspace.info[9];
      this->workspace.fact.resize(static_cast<size_t>(lfact));
      this->workspace.ifact.resize(static_cast<size_t>(lifact));
      // store the sizes of the symbolic analysis
      this->workspace.n = n;
      this->workspace.nnz = nnz;
      this->workspace.lfact = lfact;
      this->workspace.lifact = lifact;
   }

   void MA57Solver::do_numerical_factorization(const SymmetricMatrix<size_t, double>& matrix) {
      assert(matrix.dimension() <= this->dimension && "MA57Solver: the dimension of the matrix is larger than the preallocated size");
      assert(this->workspace.nnz == static_cast<int>(matrix.number_nonzeros()) && "MA57Solver: the numbers of nonzeros do not match");

      const int n = static_cast<int>(matrix.dimension());
      int nnz = static_cast<int>(matrix.number_nonzeros());

      // numerical factorization
      MA57_numerical_factorization(&n, &nnz, matrix.data_pointer(), this->workspace.fact.data(), &this->workspace.lfact,
         this->workspace.ifact.data(), &this->workspace.lifact, &this->workspace.lkeep, this->workspace.keep.data(),
         this->workspace.iwork.data(), this->workspace.icntl.data(), this->workspace.cntl.data(), this->workspace.info.data(),
         this->workspace.rinfo.data());
   }

   void MA57Solver::solve_indefinite_system(const SymmetricMatrix<size_t, double>& matrix, const Vector<double>& rhs, Vector<double>& result) {
      // solve
      const int n = static_cast<int>(matrix.dimension());
      int nnz = static_cast<int>(matrix.number_nonzeros());
      const int lrhs = n; // integer, length of rhs

      // solve the linear system
      if (this->use_iterative_refinement) {
         MA57_linear_solve_with_iterative_refinement(&this->workspace.job, &n, &nnz, matrix.data_pointer(), this->row_indices.data(),
            this->column_indices.data(), this->workspace.fact.data(), &this->workspace.lfact, this->workspace.ifact.data(),
            &this->workspace.lifact, rhs.data(), result.data(), this->workspace.residuals.data(), this->workspace.work.data(),
            this->workspace.iwork.data(), this->workspace.icntl.data(), this->workspace.cntl.data(), this->workspace.info.data(),
            this->workspace.rinfo.data());
      }
      else {
         // copy rhs into result (overwritten by MA57)
         result = rhs;

         MA57_linear_solve(&this->workspace.job, &n, this->workspace.fact.data(), &this->workspace.lfact, this->workspace.ifact.data(),
            &this->workspace.lifact, &this->workspace.nrhs, result.data(), &lrhs, this->workspace.work.data(), &this->workspace.lwork,
            this->workspace.iwork.data(), this->workspace.icntl.data(), this->workspace.info.data());
      }
   }

   void MA57Solver::solve_indefinite_system(Statistics& statistics, const Subproblem& subproblem, Direction& direction,
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
      this->solve_indefinite_system(this->augmented_matrix, this->rhs, this->solution);
      // assemble the full primal-dual direction
      subproblem.assemble_primal_dual_direction(this->solution, direction);
   }

   Inertia MA57Solver::get_inertia() const {
      // rank = number_positive_eigenvalues + number_negative_eigenvalues
      // n = rank + number_zero_eigenvalues
      const size_t rank = this->rank();
      const size_t number_negative_eigenvalues = this->number_negative_eigenvalues();
      const size_t number_positive_eigenvalues = rank - number_negative_eigenvalues;
      const size_t number_zero_eigenvalues = static_cast<size_t>(this->workspace.n) - rank;
      return {number_positive_eigenvalues, number_negative_eigenvalues, number_zero_eigenvalues};
   }

   size_t MA57Solver::number_negative_eigenvalues() const {
      return static_cast<size_t>(this->workspace.info[23]);
   }

   /*
   bool MA57Solver::matrix_is_positive_definite() const {
      // positive definite = non-singular and no negative eigenvalues
      return not this->matrix_is_singular() && this->number_negative_eigenvalues() == 0;
   }
   */

   bool MA57Solver::matrix_is_singular() const {
      return (this->workspace.info[0] == 4);
   }

   size_t MA57Solver::rank() const {
      return static_cast<size_t>(this->workspace.info[24]);
   }

   // protected member functions

   void MA57Solver::save_sparsity_pattern_internally(const SymmetricMatrix<size_t, double>& matrix) {
      // build the internal matrix representation
      this->row_indices.clear();
      this->column_indices.clear();
      for (const auto [row_index, column_index, _]: matrix) {
         this->row_indices.emplace_back(static_cast<int>(row_index + MA57Solver::fortran_shift));
         this->column_indices.emplace_back(static_cast<int>(column_index + MA57Solver::fortran_shift));
      }
   }
} // namespace