#ifndef LBFGSB_H
#define LBFGSB_H

#include "Problem.hpp"
#include "Iterate.hpp"
#include "LocalSolution.hpp"

class LBFGSB {
    public:
        LBFGSB(int memory_size = 5);
        void initialize(std::map<int,int> slacked_constraints);
        LocalSolution solve(Problem& problem, Iterate& current_point);
        
        // penalty parameter
        double rho;
        int memory_size;
    
    private:
        double compute_augmented_lagrangian_(Problem& problem, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers);
        std::vector<double> compute_augmented_lagrangian_gradient_(Problem& problem, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers);
        
        /* map of (constraint index, slack index) */
        std::map<int,int> slacked_constraints_;
        
        /* Fortran parameters */
        char task_[60];
        char csave_[60];
        int lsave_[4];
        int isave_[44];
        double dsave_[29];
        int iprint_ = 1;
        double factr_ = 1e5;
        double pgtol_ = 1e-5;
};

#endif // LBFGSB_H
