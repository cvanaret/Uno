#ifndef HESSIANEVALUATION_H
#define HESSIANEVALUATION_H

#include <ostream>
#include <map>
#include <vector>
#include "Problem.hpp"
#include "Iterate.hpp"
#include "Matrix.hpp"

// virtual (abstract) class
class HessianEvaluation {
public:
    HessianEvaluation(int size);
    virtual ~HessianEvaluation();

    int size;
    virtual void compute(Problem& problem, Iterate& iterate) = 0;
};

class ExactHessianEvaluation : public HessianEvaluation {
    /* Coordinate list */
public:
    ExactHessianEvaluation(int size);

    void compute(Problem& problem, Iterate& iterate);
};

class BFGSHessianEvaluation : public HessianEvaluation {
    /* Coordinate list */
public:
    BFGSHessianEvaluation(int size);

    void compute(Problem& problem, Iterate& iterate);
    
private:
    CSCMatrix previous_hessian;
    std::vector<double> previous_x;
};

#endif // HESSIANEVALUATION_H
