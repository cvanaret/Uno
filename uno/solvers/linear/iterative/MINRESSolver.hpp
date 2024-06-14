// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MINRESSOLVER_H
#define UNO_MINRESSOLVER_H

#include <vector>
#include <cmath>
#include "linear_algebra/Vector.hpp"
#include "solvers/linear/SymmetricIndefiniteLinearSolver.hpp"
#include "symbolic/Expression.hpp"

enum class Status {
   INITIAL_POINT_OPTIMAL,
   OPTIMAL_SOLUTION,
   LEAST_SQUARES_SOLUTION,
   REASONABLE_ACCURACY,
   EIGENVECTOR,
   LARGE_CONDITION_NUMBER,
   ITERATION_LIMIT,
   RHS_IS_EIGENVECTOR,
   USER_DEFINED
};

template <typename NumericalType>
struct GivensRotation {
   NumericalType cosine;
   NumericalType sine;
};

// implementation of MINRES based on https://web.stanford.edu/group/SOL/dissertations/sou-cheng-choi-thesis.pdf, p26
template <typename IndexType, typename NumericalType, typename LinearOperator>
class MINRESSolver : public SymmetricIndefiniteLinearSolver<IndexType, NumericalType> {
public:
   MINRESSolver(const LinearOperator& linear_operator, size_t max_dimension);
   ~MINRESSolver() override = default;

   void solve_indefinite_system(const SymmetricMatrix<IndexType, NumericalType>& matrix, const Vector<NumericalType>& rhs,
         Vector<NumericalType>& result, bool from_scratch) override;

private:
   static constexpr NumericalType tolerance{1e-6};
   const LinearOperator& linear_operator;
   const size_t max_dimension;
   size_t iteration{0};
   Vector<NumericalType> p_k;
   Vector<NumericalType> v_km1;
   Vector<NumericalType> v_k;
   Vector<NumericalType> d_km2;
   Vector<NumericalType> d_km1;
   Vector<NumericalType> d_k;
   Vector<NumericalType> x_k;
   NumericalType epsilon_k{0};
   NumericalType phi_km1;
   NumericalType delta1_k{0};
   GivensRotation<NumericalType> givens_rotation{-1, 0};
   
   // functions
   [[nodiscard]] std::pair<NumericalType, NumericalType> compute_lanczos_step(NumericalType beta_k);
   [[nodiscard]] NumericalType compute_symmetric_orthogonalization(NumericalType a, NumericalType b);
};

template <typename IndexType, typename NumericalType, typename LinearOperator>
MINRESSolver<IndexType, NumericalType, LinearOperator>::MINRESSolver(const LinearOperator& linear_operator, size_t max_dimension):
   SymmetricIndefiniteLinearSolver<IndexType, NumericalType>(max_dimension),
   linear_operator(linear_operator),
   max_dimension(max_dimension),
   p_k(max_dimension),
   v_km1(max_dimension), v_k(max_dimension),
   d_km2(max_dimension), d_km1(max_dimension), d_k(max_dimension),
   x_k(max_dimension) {
}

// TODO: user-defined termination criterion
template <typename IndexType, typename NumericalType, typename LinearOperator>
void MINRESSolver<IndexType, NumericalType, LinearOperator>::solve_indefinite_system(const SymmetricMatrix<IndexType, NumericalType>& /*matrix*/,
      const Vector<NumericalType>& rhs, Vector<NumericalType>& result, bool /*from_scratch*/) {
   this->iteration = 1;
   // v_0 = 0
   this->v_km1.fill(NumericalType(0));
   // compute beta_1 and v_1
   NumericalType beta_k = std::sqrt(rhs.dot(rhs));
   // set phi_0 = tau_0 = beta_1
   this->phi_km1 = beta_k;
   
   this->v_k = rhs;
   this->v_k.scale(NumericalType(1)/beta_k);
   
   bool termination = false;
   while (not termination) {
      // compute the Lanczos step
      const auto [alpha_k, beta_kp1] = this->compute_lanczos_step(beta_k);
      
      // last left orthogonalization on middle two entries in last column of Tk
      const NumericalType delta2_k = this->givens_rotation.cosine * this->delta1_k + this->givens_rotation.sine * alpha_k;
      const NumericalType gamma1_k = this->givens_rotation.sine * this->delta1_k - this->givens_rotation.cosine * alpha_k;
      
      // test for negative curvature descent
      if (this->givens_rotation.cosine * gamma1_k >= NumericalType(0)) {
         // std::cout << "DIRECTION OF NEGATIVE CURVATURE\n";
      }
      
      // last left orthogonalization to produce first two entries of T_k+1 e_k+1
      const NumericalType epsilon_kp1 = this->givens_rotation.sine * beta_kp1;
      const NumericalType delta1_kp1 = -this->givens_rotation.cosine * beta_kp1;
      
      // current left orthogonalization to zero out beta_k+1
      const NumericalType gamma2_k = this->compute_symmetric_orthogonalization(gamma1_k, beta_kp1);
      
      // right-hand side, residual norms, and matrix norm
      const NumericalType tau_k = this->givens_rotation.cosine * this->phi_km1; // step length
      const NumericalType phi_k = this->givens_rotation.sine * this->phi_km1;
      // TODO
      
      // update solution and matrix condition number
      if (std::abs(gamma2_k) >= this->tolerance) {
         // update d_k
         this->d_k = this->v_k - delta2_k * this->d_km1 - this->epsilon_k * this->d_km2;
         this->d_k.scale(NumericalType(1)/gamma2_k);
         
         // update x_k
         this->x_k += tau_k * this->d_k;
         
         // update v_k
         this->v_km1 = this->v_k;
         if (std::abs(beta_kp1) >= this->tolerance) {
            this->v_k = this->p_k;
            this->v_k.scale(NumericalType(1)/beta_kp1);
         }
         else {
            break;
         }
         // update the quantities for the next iteration
         beta_k = beta_kp1;
         this->delta1_k = delta1_kp1; // NumericalType
         this->epsilon_k = epsilon_kp1; // NumericalType
         this->phi_km1 = phi_k; // NumericalType
         this->d_km2 = this->d_km1; // Vector<NumericalType>
         this->d_km1 = this->d_k; // Vector<NumericalType>
      }
      else {
         break;
      }
      this->iteration++;
   }
   result = this->x_k;
}

template <typename IndexType, typename NumericalType, typename LinearOperator>
std::pair<NumericalType, NumericalType> MINRESSolver<IndexType, NumericalType, LinearOperator>::compute_lanczos_step(NumericalType beta_k) {
   // compute matrix-vector product in p_k
   this->linear_operator(this->v_k, this->p_k);
   
   // compute alpha_k = v_k^T A v_k
   const NumericalType alpha_k = this->v_k.dot(this->p_k);
   
   // compute p_k = p_k - beta_k v_k-1 - alpha_k v_k
   this->p_k = this->p_k - beta_k * this->v_km1 - alpha_k * this->v_k;
   
   // compute beta_k+1 = ||p_k||
   const NumericalType beta_kp1 = std::sqrt(this->p_k.dot(this->p_k));
   return {alpha_k, beta_kp1};
}

template <typename NumericalType>
NumericalType sign(NumericalType x) {
   if (x < NumericalType(0)) {
      return NumericalType(-1);
   }
   else if (x > NumericalType(0)) {
      return NumericalType(1);
   }
   else {
      return NumericalType(0);
   }
}

template <typename IndexType, typename NumericalType, typename LinearOperator>
NumericalType MINRESSolver<IndexType, NumericalType, LinearOperator>::compute_symmetric_orthogonalization(NumericalType a, NumericalType b) {
   if (b == NumericalType(0)) {
      this->givens_rotation.sine = NumericalType(0);
      if (a == NumericalType(0)) {
         this->givens_rotation.cosine = NumericalType(1);
      }
      else {
         this->givens_rotation.cosine = sign(a);
      }
      return std::abs(a);
   }
   else if (a == NumericalType(0)) {
      this->givens_rotation = {NumericalType(0), sign(b)};
      return std::abs(b);
   }
   else if (std::abs(b) > std::abs(a)) {
      const NumericalType tau = a/b;
      const NumericalType sine = sign(b)/std::sqrt(1 + tau*tau);
      const NumericalType cosine = sine*tau;
      this->givens_rotation = {cosine, sine};
      return b/sine;
   }
   else {
      const NumericalType tau = b/a;
      const NumericalType cosine = sign(a)/std::sqrt(1 + tau*tau);
      const NumericalType sine = cosine*tau;
      this->givens_rotation = {cosine, sine};
      return a/cosine;
   }
}

#endif // UNO_MINRESSOLVER_H
