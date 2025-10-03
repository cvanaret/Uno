// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_AMPLMODEL_H
#define UNO_AMPLMODEL_H

#include <vector>
#include "model/Model.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"
#include "symbolic/CollectionAdapter.hpp"
// include AMPL Solver Library (ASL)
extern "C" {
#include "asl_pfgh.h"
#include "getstub.h"
}

namespace uno {
   // forward declaration
   class Result;

   class AMPLModel: public Model {
   public:
      explicit AMPLModel(const std::string& file_name);

      static constexpr double lagrangian_sign_convention{-1.};

      ~AMPLModel() override;

      // availability of linear operators
      [[nodiscard]] bool has_jacobian_operator() const override;
      [[nodiscard]] bool has_jacobian_transposed_operator() const override;
      [[nodiscard]] bool has_hessian_operator() const override;
      [[nodiscard]] bool has_hessian_matrix() const override;

      // function evaluations
      [[nodiscard]] double evaluate_objective(const Vector<double>& x) const override;
      void evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const override;

      // dense objective gradient
      void evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const override;

      // structures of Jacobian and Hessian
      void compute_constraint_jacobian_sparsity(int* row_indices, int* column_indices, int solver_indexing,
         MatrixOrder matrix_order) const override;
      void compute_hessian_sparsity(int* row_indices, int* column_indices, int solver_indexing) const override;

      // numerical evaluations of Jacobian and Hessian
      void evaluate_constraint_jacobian(const Vector<double>& x, double* jacobian_values) const override;
      void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         double* hessian_values) const override;

      // linear operators for Jacobian-, Jacobian^T-, and Hessian-vector products
      void compute_jacobian_vector_product(const double* x, const double* vector, double* result) const override;
      void compute_jacobian_transposed_vector_product(const double* x, const double* vector, double* result) const override;
      void compute_hessian_vector_product(const double* x, const double* vector, double objective_multiplier,
         const Vector<double>& multipliers, double* result) const override;

      [[nodiscard]] double variable_lower_bound(size_t variable_index) const override;
      [[nodiscard]] double variable_upper_bound(size_t variable_index) const override;
      [[nodiscard]] const SparseVector<size_t>& get_slacks() const override;
      [[nodiscard]] const Vector<size_t>& get_fixed_variables() const override;

      [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override;
      [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override;
      [[nodiscard]] const Collection<size_t>& get_equality_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_inequality_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_linear_constraints() const override;

      void initial_primal_point(Vector<double>& x) const override;
      void initial_dual_point(Vector<double>& multipliers) const override;
      void postprocess_solution(Iterate& iterate) const override;

      void write_solution_to_file(Result& result) const;

      [[nodiscard]] size_t number_jacobian_nonzeros() const override;
      [[nodiscard]] size_t number_hessian_nonzeros() const override;

   private:
      // private constructor to pass the dimensions to the Model base constructor
      AMPLModel(const std::string& file_name, ASL* asl);

      // mutable: can be modified by const methods (internal state not seen by user)
      mutable ASL* asl; /*!< Instance of the AMPL Solver Library class */
      size_t number_asl_hessian_nonzeros{0}; /*!< Number of nonzero elements in the Hessian */

      // lists of variables and constraints + corresponding collection objects
      ForwardRange linear_constraints;
      std::vector<size_t> equality_constraints{};
      CollectionAdapter<std::vector<size_t>&> equality_constraints_collection;
      std::vector<size_t> inequality_constraints{};
      CollectionAdapter<std::vector<size_t>&> inequality_constraints_collection;
      SparseVector<size_t> slacks{};
      Vector<size_t> fixed_variables;

      void compute_lagrangian_hessian_sparsity();
   };

   // check that an array of integers is in increasing order (x[i] <= x[i+1])
   template <typename Array>
   bool in_increasing_order(const Array& array, size_t length) {
      size_t index = 0;
      while (index < length - 1) {
         if (array[index] > array[index + 1]) {
            return false;
         }
         index++;
      }
      return true;
   }
} // namespace

#endif // UNO_AMPLMODEL_H