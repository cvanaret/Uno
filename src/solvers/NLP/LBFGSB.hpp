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
        LocalSolution solve(Problem& problem, Iterate& current_iterate,
            double (*compute_objective)(Problem&, std::vector<double>&, std::vector<double>&, std::vector<double>&, double),
            std::vector<double> (*compute_objective_gradient)(Problem&, std::map<int,int>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, double),
            std::vector<double> (*compute_constraints)(Problem& problem, std::map<int,int>& slacked_constraints, std::vector<double>& x),
            std::vector<double>& l, std::vector<double>& u, std::vector<ConstraintType>& variable_status,
            int max_iterations);
        
        // TODO REMOVE augmented Lagrangian penalty parameter
        double penalty_parameter;
        int limited_memory_size;  // number limited memory vectors (= m in lbfgsb.f)
    
    private:
        /* TODO REMOVE map of (constraint index, slack index) */
        std::map<int,int> slacked_constraints_;
        
        /* Fortran parameters needed by lbfgsb.f */
        char task_[60];
        char csave_[60];
        int lsave_[4];
        int isave_[44];
        double dsave_[29];
        int iprint_ = -1;
        double factr_ = 1e5;
        double pgtol_ = 1e-5;
};

#endif // LBFGSB_H
