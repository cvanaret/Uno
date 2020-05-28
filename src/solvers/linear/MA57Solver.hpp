#ifndef MA57SOLVER_H
#define MA57SOLVER_H

#include <vector>
#include "Matrix.hpp"
#include "LinearSolver.hpp"

struct MA57Factorization {
    int dimension;
    std::vector<double> fact;
    int lfact;
    std::vector<int> ifact;
    int lifact;
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

    short use_fortran;

    void factorize(COOMatrix& matrix) override;
    void solve(std::vector<double>& rhs) override;
    int number_negative_eigenvalues() override;
    bool matrix_is_singular() override;
    int rank() override;

private:
    /* for ma57id_ */
    std::vector<double> cntl_;
    std::vector<int> icntl_;
    std::vector<double> rinfo_;
    
    MA57Factorization factorization_;
};

#endif // MA57SOLVER_H
