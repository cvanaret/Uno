#include "ipopt_uno_adapter.h"
#include "../Uno_C_API.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// --------------------
// Internal problem struct
// --------------------
struct IpoptProblemInfo {
    ipindex n;           // Number of variables
    ipindex m;           // Number of constraints
    ipindex nele_jac;
    ipindex nele_hess;
    ipindex index_style;

    // Bounds
    ipnumber* x_L;
    ipnumber* x_U;
    ipnumber* g_L;
    ipnumber* g_U;

    // Callbacks
    Eval_F_CB eval_f;
    Eval_Grad_F_CB eval_grad_f;
    Eval_G_CB eval_g;
    Eval_Jac_G_CB eval_jac_g;
    Eval_H_CB eval_h;
    Intermediate_CB intermediate_cb;

    UserDataPtr user_data;

    void* uno_model;
    void* uno_solver;
};

// --------------------
// Callback wrappers
// --------------------
static uno_int objective_wrapper(uno_int n, const double* x, double* obj_val, void* user_data) {
    IpoptProblemInfo* prob = (IpoptProblemInfo*)user_data;
    return prob->eval_f(n, (ipnumber*)x, true, obj_val, prob->user_data) ? 0 : 1;
}

static uno_int gradient_wrapper(uno_int n, const double* x, double* grad, void* user_data) {
    IpoptProblemInfo* prob = (IpoptProblemInfo*)user_data;
    return prob->eval_grad_f(n, (ipnumber*)x, true, grad, prob->user_data) ? 0 : 1;
}

static uno_int constraints_wrapper(uno_int n, uno_int m, const double* x, double* g, void* user_data) {
    IpoptProblemInfo* prob = (IpoptProblemInfo*)user_data;
    return prob->eval_g(n, (ipnumber*)x, true, m, (ipnumber*)g, prob->user_data) ? 0 : 1;
}

static uno_int jacobian_wrapper(uno_int n, uno_int nnz, const double* x, double* jac, void* user_data) {
    IpoptProblemInfo* prob = (IpoptProblemInfo*)user_data;
    return prob->eval_jac_g(n, (ipnumber*)x, true, prob->m, nnz, NULL, NULL, jac, prob->user_data) ? 0 : 1;
}

static uno_int hessian_wrapper(uno_int n, uno_int m, uno_int nnz,
                               const double* x, double obj_factor, const double* lambda,
                               double* hess, void* user_data) {
    IpoptProblemInfo* prob = (IpoptProblemInfo*)user_data;
    return prob->eval_h(n, (ipnumber*)x, true, obj_factor, m, (ipnumber*)lambda,
                        true, nnz, NULL, NULL, hess, prob->user_data) ? 0 : 1;
}

// --------------------
// Internal free helper
// --------------------
static void free_prob_struct(IpoptProblemInfo* prob) {
    if(!prob) return;
    free(prob->x_L); free(prob->x_U);
    free(prob->g_L); free(prob->g_U);
    if(prob->uno_solver) uno_destroy_solver(prob->uno_solver);
    if(prob->uno_model) uno_destroy_model(prob->uno_model);
    free(prob);
}

// --------------------
// Create / Free
// --------------------
IpoptProblem CreateIpoptProblem(ipindex n, const ipnumber* x_L, const ipnumber* x_U,
                                ipindex m, const ipnumber* g_L, const ipnumber* g_U,
                                ipindex nele_jac, ipindex nele_hess, ipindex index_style,
                                Eval_F_CB eval_f, Eval_G_CB eval_g, Eval_Grad_F_CB eval_grad_f,
                                Eval_Jac_G_CB eval_jac_g, Eval_H_CB eval_h)
{
    IpoptProblemInfo* prob = (IpoptProblemInfo*)malloc(sizeof(IpoptProblemInfo));
    if(!prob) return NULL;

    prob->n = n; prob->m = m; prob->nele_jac = nele_jac; prob->nele_hess = nele_hess;
    prob->index_style = index_style;

    prob->x_L = (ipnumber*)malloc(sizeof(ipnumber)*n);
    prob->x_U = (ipnumber*)malloc(sizeof(ipnumber)*n);
    prob->g_L = (ipnumber*)malloc(sizeof(ipnumber)*m);
    prob->g_U = (ipnumber*)malloc(sizeof(ipnumber)*m);

    if(!prob->x_L || !prob->x_U || !prob->g_L || !prob->g_U) {
        free_prob_struct(prob);
        return NULL;
    }

    memcpy(prob->x_L, x_L, sizeof(ipnumber)*n);
    memcpy(prob->x_U, x_U, sizeof(ipnumber)*n);
    memcpy(prob->g_L, g_L, sizeof(ipnumber)*m);
    memcpy(prob->g_U, g_U, sizeof(ipnumber)*m);

    prob->eval_f = eval_f; prob->eval_grad_f = eval_grad_f;
    prob->eval_g = eval_g; prob->eval_jac_g = eval_jac_g;
    prob->eval_h = eval_h; prob->intermediate_cb = NULL;
    prob->user_data = NULL;

    // create Uno model
    prob->uno_model = uno_create_model(UNO_PROBLEM_NONLINEAR, n, prob->x_L, prob->x_U, index_style);
    if(!prob->uno_model) {
        free_prob_struct(prob);
        return NULL;
    }
    uno_set_objective(prob->uno_model, UNO_MINIMIZE, objective_wrapper, gradient_wrapper);
    if(m>0)
        uno_set_constraints(prob->uno_model, m, constraints_wrapper, prob->g_L, prob->g_U, nele_jac, NULL, NULL, jacobian_wrapper);

    // create Uno solver
    prob->uno_solver = uno_create_solver();
    if(!prob->uno_solver) {
        free_prob_struct(prob);
        return NULL;
    }

    return prob;
}

void FreeIpoptProblem(IpoptProblem prob)
{
    free_prob_struct(prob);
}

// --------------------
// Options
// --------------------
bool AddIpoptStrOption(IpoptProblem prob, const char* keyword, const char* val)
{
    return uno_set_solver_string_option(prob->uno_solver, keyword, val);
}

bool AddIpoptNumOption(IpoptProblem prob, const char* keyword, ipnumber val)
{
    return uno_set_solver_double_option(prob->uno_solver, keyword, val);
}

bool AddIpoptIntOption(IpoptProblem prob, const char* keyword, ipindex val)
{
    return uno_set_solver_integer_option(prob->uno_solver, keyword, val);
}

bool SetIntermediateCallback(IpoptProblem prob, Intermediate_CB cb)
{
    prob->intermediate_cb = cb;
    return true;
}

// --------------------
// Solve
// --------------------
ApplicationReturnStatus IpoptSolve(IpoptProblem prob, ipnumber* x, ipnumber* g,
                                   ipnumber* obj_val, ipnumber* mult_g,
                                   ipnumber* mult_x_L, ipnumber* mult_x_U,
                                   UserDataPtr user_data)
{
    prob->user_data = user_data;
    if(!prob->uno_solver || !prob->uno_model) return Unknown_Error;

    uno_set_user_data(prob->uno_model, prob);
    uno_optimize(prob->uno_solver, prob->uno_model);

    for(ipindex i=0;i<prob->n;i++)
        x[i] = uno_get_primal_solution_component(prob->uno_solver, i);
    if(g) for(ipindex i=0;i<prob->m;i++)
        g[i] = 0.0; // only constraint evaluation placeholder

    *obj_val = uno_get_solution_objective(prob->uno_solver);

    for(ipindex i=0;i<prob->n;i++){
        if(mult_x_L) mult_x_L[i] = uno_get_lower_bound_dual_solution_component(prob->uno_solver, i);
        if(mult_x_U) mult_x_U[i] = uno_get_upper_bound_dual_solution_component(prob->uno_solver, i);
    }
    if(mult_g) for(ipindex i=0;i<prob->m;i++)
        mult_g[i] = uno_get_constraint_dual_solution_component(prob->uno_solver, i);

    uno_int status = uno_get_optimization_status(prob->uno_solver);
    switch(status){
        case UNO_SUCCESS: return Solve_Succeeded;
        case UNO_USER_TERMINATION: return User_Requested_Stop;
        case UNO_ITERATION_LIMIT: return Maximum_Iterations_Exceeded;
        default: return Unknown_Error;
    }
}

// --------------------
// Version
// --------------------
void GetIpoptVersion(int* major, int* minor, int* release)
{
    if(!major || !minor || !release) return;
    uno_get_version((uno_int*)major, (uno_int*)minor, (uno_int*)release);
}
