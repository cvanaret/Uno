// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MINRESSOLVER_H
#define UNO_MINRESSOLVER_H

#include <vector>
#include "linear_algebra/Vector.hpp"
#include "solvers/linear/SymmetricIndefiniteLinearSolver.hpp"

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

template <typename NumericalType, typename LinearOperator>
class MINRESSolver : public SymmetricIndefiniteLinearSolver<NumericalType> {
public:
   explicit MINRESSolver(const LinearOperator& linear_operator, size_t max_dimension);
   ~MINRESSolver() override = default;

   void factorize(const SymmetricMatrix<NumericalType>& /*matrix*/) override { /* do nothing */  }
   void do_symbolic_factorization(const SymmetricMatrix<NumericalType>& /*matrix*/) override { /* do nothing */ }
   void do_numerical_factorization(const SymmetricMatrix<NumericalType>& /*matrix*/) override { /* do nothing */ }
   void solve_indefinite_system(const SymmetricMatrix<NumericalType>& matrix, const std::vector<NumericalType>& rhs, std::vector<NumericalType>& result) override;

   [[nodiscard]] std::tuple<size_t, size_t, size_t> get_inertia() const override { return {0, 0, 0}; }
   [[nodiscard]] size_t number_negative_eigenvalues() const override { return 0; }
   [[nodiscard]] bool matrix_is_singular() const override { return false; }
   [[nodiscard]] size_t rank() const override { return 0; }

private:
   static constexpr NumericalType tolerance{1e-6};
   const LinearOperator& linear_operator;
   const size_t max_dimension;
   size_t iteration{0};
   std::vector<NumericalType> p_k;
   std::vector<NumericalType> v_km1;
   std::vector<NumericalType> v_k;
   std::vector<NumericalType> v_kp1;
   std::vector<NumericalType> d_km2;
   std::vector<NumericalType> d_km1;
   std::vector<NumericalType> d_k;
   std::vector<NumericalType> x_k;
   //NumericalType beta_k;
   NumericalType epsilon_k{0};
   NumericalType phi_km1;
   NumericalType tau_km1;
   NumericalType delta1_k{0};
   GivensRotation<NumericalType> givens_rotation_km1{-1, 0};
   
   // functions
   std::pair<NumericalType, NumericalType> compute_lanczos_step(const LinearOperator& linear_operator, NumericalType beta_k);
   NumericalType compute_symmetric_orthogonalization(NumericalType a, NumericalType b);
};

template <typename NumericalType, typename LinearOperator>
MINRESSolver<NumericalType, LinearOperator>::MINRESSolver(const LinearOperator& linear_operator, size_t max_dimension):
   SymmetricIndefiniteLinearSolver<NumericalType>(max_dimension),
   linear_operator(linear_operator),
   max_dimension(max_dimension),
   p_k(max_dimension),
   v_km1(max_dimension), v_k(max_dimension), v_kp1(max_dimension),
   d_km2(max_dimension), d_km1(max_dimension), d_k(max_dimension),
   x_k(max_dimension) {
}

// TODO: user-defined termination criterion
template <typename NumericalType, typename LinearOperator>
void MINRESSolver<NumericalType, LinearOperator>::solve_indefinite_system(const SymmetricMatrix<NumericalType>& /*matrix*/,
      const std::vector<NumericalType>& rhs, std::vector<NumericalType>& result) {
   this->iteration = 1;
   // v_0 = 0
   initialize_vector(this->v_km1, NumericalType(0));
   // compute beta_1 and v_1
   NumericalType beta_k = std::sqrt(dot_product(rhs, rhs));
   // set phi_0 = tau_0 = beta_1
   this->phi_km1 = this->tau_km1 = beta_k;
   // std::cout << "beta_1 = " << beta_k << '\n';
   
   this->v_k = rhs;
   scale(this->v_k, NumericalType(1)/beta_k);
   
   bool termination = false;
   while (not termination) {
      // std::cout << "\nCurrent x: "; print_vector(std::cout, this->x_k);
      // std::cout << "v_" << this->iteration << " = "; print_vector(std::cout, v_k);
      
      // compute the Lanczos step
      const auto [alpha_k, beta_kp1] = this->compute_lanczos_step(this->linear_operator, beta_k);
      // std::cout << "beta_" << this->iteration+1 << " = " << lanczos_coefs.beta_kp1 << '\n';
      
      // last left orthogonalization on middle two entries in last column of Tk
      const NumericalType delta2_k = this->givens_rotation_km1.cosine*this->delta1_k + this->givens_rotation_km1.sine*alpha_k;
      const NumericalType gamma1_k = this->givens_rotation_km1.sine*this->delta1_k - this->givens_rotation_km1.cosine*alpha_k;
      
      // test for negative curvature descent
      if (this->givens_rotation_km1.cosine*gamma1_k >= NumericalType(0)) {
         // std::cout << this->givens_rotation_km1.cosine*gamma1_k << '\n';
         // std::cout << "DIRECTION OF NEGATIVE CURVATURE\n";
      }
      
      // last left orthogonalization to produce first two entries of T_k+1 e_k+1
      const NumericalType epsilon_kp1 = this->givens_rotation_km1.sine*beta_kp1;
      const NumericalType delta1_kp1 = -this->givens_rotation_km1.cosine*beta_kp1;
      
      // current left orthogonalization to zero out beta_k+1
      const NumericalType gamma2_k = this->compute_symmetric_orthogonalization(gamma1_k, beta_kp1);
      
      // right-hand side, residual norms, and matrix norm
      const NumericalType tau_k = this->givens_rotation_km1.cosine*this->phi_km1; // step length
      const NumericalType phi_k = this->givens_rotation_km1.sine*this->phi_km1;
      // TODO
      
      // update solution and matrix condition number
      if (std::abs(gamma2_k) >= NumericalType(tolerance)) {
         // update d_k
         add_vectors(NumericalType(1), this->v_k, -delta2_k, this->d_km1, this->d_k);
         add_vectors(NumericalType(1), this->d_k, -this->epsilon_k, this->d_km2, this->d_k);
         scale(this->d_k, NumericalType(1)/gamma2_k);
         
         // update x_k
         add_vectors(NumericalType(1), this->x_k, tau_k, this->d_k, this->x_k);
         
         // update v_k
         if (std::abs(beta_kp1) >= NumericalType(tolerance)) {
            this->v_kp1 = this->p_k;
            scale(this->v_kp1, NumericalType(1)/beta_kp1);
         }
         else {
            break;
         }
         beta_k = beta_kp1;
         this->delta1_k = delta1_kp1;
         this->epsilon_k = epsilon_kp1;
         this->tau_km1 = tau_k;
         this->phi_km1 = phi_k;
         this->tau_km1 = tau_k;
         
         this->v_km1 = this->v_k;
         this->v_k = this->v_kp1;
         this->d_km2 = this->d_km1;
         this->d_km1 = this->d_k;
      }
      else {
         break;
      }
      this->iteration++;
   }
   copy_from(result, this->x_k);
}

template <typename NumericalType, typename LinearOperator>
std::pair<NumericalType, NumericalType> MINRESSolver<NumericalType, LinearOperator>::compute_lanczos_step(const LinearOperator& linear_operator,
      NumericalType beta_k) {
   // std::cout << "Starting compute_lanczos_step\n";
   // compute matrix-vector product in p_k
   linear_operator(this->v_k, this->p_k);
   
   // compute alpha_k = v_k^T A v_k
   const NumericalType alpha_k = dot_product(this->v_k, this->p_k);
   
   // compute the full p_k = p_k - beta_k v_k-1 - alpha_k v_k
   add_vectors(NumericalType(1), this->p_k, -alpha_k, this->v_k, this->p_k);
   add_vectors(NumericalType(1), this->p_k, -beta_k, this->v_km1, this->p_k);
   
   // compute beta_k+1 = ||p_k||
   const NumericalType beta_kp1 = std::sqrt(dot_product(this->p_k, this->p_k));
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

template <typename NumericalType, typename LinearOperator>
NumericalType MINRESSolver<NumericalType, LinearOperator>::compute_symmetric_orthogonalization(NumericalType a, NumericalType b) {
   if (b == NumericalType(0)) {
      this->givens_rotation_km1.sine = NumericalType(0);
      if (a == NumericalType(0)) {
         this->givens_rotation_km1.cosine = NumericalType(1);
      }
      else {
         this->givens_rotation_km1.cosine = sign(a);
      }
      return std::abs(a);
   }
   else if (a == NumericalType(0)) {
      this->givens_rotation_km1 = {NumericalType(0), sign(b)};
      return std::abs(b);
   }
   else if (std::abs(b) > std::abs(a)) {
      const NumericalType tau = a/b;
      const NumericalType sine = sign(b)/std::sqrt(1 + tau*tau);
      const NumericalType cosine = sine*tau;
      this->givens_rotation_km1 = {cosine, sine};
      return b/sine;
   }
   else {
      const NumericalType tau = b/a;
      const NumericalType cosine = sign(a)/std::sqrt(1 + tau*tau);
      const NumericalType sine = cosine*tau;
      this->givens_rotation_km1 = {cosine, sine};
      return a/cosine;
   }
}

#endif // UNO_MINRESSOLVER_H
