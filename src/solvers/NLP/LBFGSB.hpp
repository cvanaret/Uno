#ifndef LBFGSB_H
#define LBFGSB_H

#include "Problem.hpp"
#include "Iterate.hpp"
#include "LocalSolution.hpp"

/*! \class LBFGSB
 * \brief Limited-memory BFGS for Bounds (L-BFGS-B)
 *
 *  interface to bound-constrained nonlinear optimization solver L-BFGS-B by Nocedal et al.
 */
class LBFGSB {
    public:
        LBFGSB(int limited_memory_size = 5);
        void initialize(std::map<int,int> slacked_constraints);
        LocalSolution solve(Problem& problem, Iterate& current_point);
        LocalSolution solve(Problem& problem, Iterate& current_iterate,
            double (*compute_augmented_lagrangian)(Problem&, std::map<int,int>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, double),
            std::vector<double> (*compute_augmented_lagrangian_gradient)(Problem&, std::map<int,int>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, double));
        
        // augmented Lagrangian penalty parameter
        double rho;
        int limited_memory_size;  // number limited memory vectors (=m in lbfgsb.f)
    
    private:
        // TODO move to AL
        double compute_augmented_lagrangian_(Problem& problem, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers);
        std::vector<double> compute_augmented_lagrangian_gradient_(Problem& problem, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers);
        /* map of (constraint index, slack index) */
        std::map<int,int> slacked_constraints_;
        
        /* Fortran parameters needed by lbfgsb.f */
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
