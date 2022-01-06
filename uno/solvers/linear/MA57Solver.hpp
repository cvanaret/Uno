#ifndef UNO_MA57SOLVER_H
#define UNO_MA57SOLVER_H

#include <vector>
#include "LinearSolver.hpp"
#include "linear_algebra/COOSymmetricMatrix.hpp"

struct MA57Factorization {
   int n{};
   int nnz{};
   std::vector<double> fact{};
   int lfact{};
   std::vector<int> ifact{};
   int lifact{};
   int lkeep{};
   std::vector<int> keep{};

   MA57Factorization() = default;
};

/*! \class MA57Solver
 * \brief Interface for MA57
 * see https://github.com/YimingYAN/linSolve
 *
 *  Interface to the sparse symmetric linear solver MA57
 */
class MA57Solver : public LinearSolver {
public:
   explicit MA57Solver(size_t max_dimension, size_t max_number_nonzeros);
   ~MA57Solver() override = default;

   void factorize(const SymmetricMatrix& matrix) override;
   void do_symbolic_factorization(const SymmetricMatrix& matrix) override;
   void do_numerical_factorization(const SymmetricMatrix& matrix) override;
   void solve(const SymmetricMatrix& matrix, const std::vector<double>& rhs, std::vector<double>& result) override;

   [[nodiscard]] std::tuple<size_t, size_t, size_t> get_inertia() const override;
   [[nodiscard]] size_t number_negative_eigenvalues() const override;
   [[nodiscard]] bool matrix_is_positive_definite() const override;
   [[nodiscard]] bool matrix_is_singular() const override;
   [[nodiscard]] size_t rank() const override;

private:
   // internal matrix representation
   std::vector<int> row_indices;
   std::vector<int> column_indices;

   std::vector<int> iwork;
   int lwork;
   std::vector<double> work;
   /* for ma57id_ (default values of controlling parameters) */
   std::array<double, 5> cntl{};
   std::array<int, 20> icntl{};
   std::array<double, 20> rinfo{};
   std::array<int, 40> info{};
   const int nrhs{1}; // number of right hand side being solved
   const int job{1};
   std::vector<double> residuals;
   const size_t fortran_shift{1};

   MA57Factorization factorization{};
   bool use_iterative_refinement{false};
   void save_matrix_to_local_format(const SymmetricMatrix& matrix);
};

#endif // UNO_MA57SOLVER_H
