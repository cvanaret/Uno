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
   std::vector<int> info;
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
   virtual ~MA57Solver() = default;

   short use_fortran;

   void factorize(COOMatrix& matrix) override;
   void do_symbolic_factorization(const COOMatrix& matrix) override;
   void do_numerical_factorization(const COOMatrix& matrix) override;
   void solve(std::vector<double>& rhs) override;

   int number_negative_eigenvalues() const override;
   bool matrix_is_singular() const override;
   int rank() const override;

private:
   /* for ma57id_ (default values of controlling parameters) */
   std::vector<double> cntl_;
   std::vector<int> icntl_;
   std::vector<double> rinfo_;

   MA57Factorization factorization_;
};

#endif // MA57SOLVER_H
