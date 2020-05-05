#ifndef HESSIANEVALUATION_H
#define HESSIANEVALUATION_H

#include <ostream>
#include <map>
#include <memory>
#include <vector>
#include "Problem.hpp"
#include "Iterate.hpp"
#include "Matrix.hpp"

// virtual (abstract) class
class HessianEvaluation {
public:
    HessianEvaluation(int dimension);
    virtual ~HessianEvaluation();

    int dimension;
    bool convexify;
    
    virtual void compute(Problem& problem, Iterate& iterate, double objective_multiplier, std::vector<double>& constraint_multipliers) = 0;
    CSCMatrix modify_inertia(CSCMatrix& hessian);
};

class ExactHessianEvaluation : public HessianEvaluation {
    /* Coordinate list */
public:
    ExactHessianEvaluation(int dimension);

    void compute(Problem& problem, Iterate& iterate, double objective_multiplier, std::vector<double>& constraint_multipliers);
};

class BFGSHessianEvaluation : public HessianEvaluation {
    /* Coordinate list */
public:
    BFGSHessianEvaluation(int dimension);

    void compute(Problem& problem, Iterate& iterate, double objective_multiplier, std::vector<double>& constraint_multipliers);
    
private:
    CSCMatrix previous_hessian;
    std::vector<double> previous_x;
};

class HessianEvaluationFactory {
    public:
		static std::shared_ptr<HessianEvaluation> create(std::string hessian_evaluation_method, int dimension);
};
#endif // HESSIANEVALUATION_H
