#ifndef LBFGSB_H
#define LBFGSB_H

#include "Problem.hpp"
#include "Iterate.hpp"
#include "LocalSolution.hpp"

class LBFGSB {
    public:
        LocalSolution solve(Problem& problem, Iterate& current_point);

};

#endif // LBFGSB_H
