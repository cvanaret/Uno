#ifndef UNO_MA27SOLVER_H
#define UNO_MA27SOLVER_H

#include <array>
#include <vector>
#include "solvers/DirectSymmetricIndefiniteLinearSolver.hpp"

namespace uno {
// forward declaration
template <typename ElementType>
class Vector;


class MA27Solver 
    : public DirectSymmetricIndefiniteLinearSolver<size_t, double>
{
public:
   explicit MA27Solver(size_t max_dimension, size_t max_number_nonzeros);
   ~MA27Solver() override = default;

    void factorize(const SymmetricMatrix<size_t, double>& matrix) override;
    void do_symbolic_factorization(const SymmetricMatrix<size_t, double>& matrix) override;
    void do_numerical_factorization(const SymmetricMatrix<size_t, double>& matrix) override;
    void solve_indefinite_system(const SymmetricMatrix<size_t, double>& matrix, const Vector<double>& rhs, Vector<double>& result) override;

    
    [[nodiscard]] std::tuple<size_t, size_t, size_t> get_inertia() const override;
    [[nodiscard]] size_t number_negative_eigenvalues() const override;
    // [[nodiscard]] bool matrix_is_positive_definite() const override;
    [[nodiscard]] bool matrix_is_singular() const override;
    [[nodiscard]] size_t rank() const override;

private:
   int nz_max{};                    // maximal number of nonzeros entries
   int n{};			                // dimension of current factorisation (maximal value here is <= max_dimension)
   int nnz{};			            // number of nonzeros of current factorisation
   std::array<int,30> icntl{};      // integer array of length 30; integer control values
   std::array<double,5> cntl{};     // double array of length 5; double control values

   std::vector<int> irn{};          // row index of input
   std::vector<int> icn{};          // col index of input

   std::vector<int> iw{};           // integer workspace of length liw
   std::vector<int> ikeep{};        // integer array of 3*n; pivot sequence
   std::vector<int> iw1{};          // integer workspace array of length n
   int nsteps{};                    // integer, to be set by ma27
   int iflag{};                     // integer; 0 if pivot order chosen automatically; 1 if the pivot order set by ikeep
   std::array<int,20> info{};       // integer array of length 20
   double ops{};                    // double, operations count

   std::vector<double> factor{};    // data array of length la;
   int maxfrt{};                    // integer, to be set by ma27


   // bool use_iterative_refinement{false}; // Not sure how to do this with ma27
   void save_matrix_to_local_format(const SymmetricMatrix<size_t, double>& matrix);
   void repeat_factorization_after_resizing();
};

} // namespace uno
#endif // UNO_MA27SOLVER_H