#ifndef LBFGSB_H
#define LBFGSB_H

#include "Problem.hpp"
#include "Iterate.hpp"
#include "LocalSolution.hpp"

class LBFGSB {
    public:
        LBFGSB();
        void initialize(std::map<int,int> slacked_constraints);
        LocalSolution solve(Problem& problem, Iterate& current_point);
        double compute_augmented_lagrangian(Problem& problem, Iterate& current_iterate, std::vector<double>& x);
        std::vector<double> compute_augmented_lagrangian_gradient(Problem& problem, Iterate& current_iterate, std::vector<double>& x);
        
        std::map<int,int> slacked_constraints;
        // penalty parameter
        double rho;
};

#endif // LBFGSB_H
