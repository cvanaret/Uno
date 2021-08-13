#ifndef MA57SOLVER_H
#define MA57SOLVER_H

#include <vector>
#include <cassert>
#include "COOSymmetricMatrix.hpp"
#include "LinearSolver.hpp"

struct MA57Factorization {
   int n;
   int nnz;
   std::vector<double> fact;
   int lfact;
   std::vector<int> ifact;
   int lifact;
   int lkeep;
   std::vector<int> keep;
   std::vector<int> iwork;
   std::array<int, 40> info;
};

/*! \class MA57Solver
 * \brief Interface for MA57
 * see https://github.com/YimingYAN/linSolve
 *
 *  Interface to the sparse symmetric linear solver MA57
 */
class MA57Solver : public LinearSolver<COOSymmetricMatrix> {
public:
   MA57Solver();
   ~MA57Solver() override = default;

   void factorize(COOSymmetricMatrix& matrix) override;
   void do_symbolic_factorization(COOSymmetricMatrix& matrix) override;
   void do_numerical_factorization(COOSymmetricMatrix& matrix) override;
   std::vector<double> solve(COOSymmetricMatrix& matrix, const std::vector<double>& rhs) override;

   [[nodiscard]] std::tuple<int, int, int> get_inertia() const override;
   [[nodiscard]] size_t number_negative_eigenvalues() const override;
   [[nodiscard]] bool matrix_is_singular() const override;
   [[nodiscard]] int rank() const override;

private:
   /* for ma57id_ (default values of controlling parameters) */
   std::array<double, 5> cntl{};
   std::array<int, 20> icntl{};
   std::array<double, 20> rinfo{};
   const int nrhs{1}; // number of right hand side being solved
   const int job{1};

   MA57Factorization factorization;
   bool use_iterative_refinement{false};
};

#endif // MA57SOLVER_H
