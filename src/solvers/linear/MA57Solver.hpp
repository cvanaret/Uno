#ifndef MA57SOLVER_H
#define MA57SOLVER_H

#include <vector>
#include "Matrix.hpp"
#include "LinearSolver.hpp"

struct MA57Factorization {
   size_t n;
   size_t nnz;
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
class MA57Solver : public LinearSolver {
public:
   MA57Solver();
   ~MA57Solver() override = default;

   short use_fortran;

   void factorize(const COOMatrix& matrix) override;
   void do_symbolic_factorization(const COOMatrix& matrix) override;
   void do_numerical_factorization(const COOMatrix& matrix) override;
   std::vector<double> solve(const COOMatrix& matrix, std::vector<double>& rhs) override;

   [[nodiscard]] std::tuple<int, int, int> get_inertia() const override;
   [[nodiscard]] size_t number_negative_eigenvalues() const override;
   [[nodiscard]] bool matrix_is_singular() const override;
   [[nodiscard]] int rank() const override;

private:
   /* for ma57id_ (default values of controlling parameters) */
   std::vector<double> cntl_;
   std::vector<int> icntl_;
   std::vector<double> rinfo_;

   MA57Factorization factorization_;
   bool use_iterative_refinement{true};
};

#endif // MA57SOLVER_H
