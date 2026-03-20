#include "ipopt_uno_adapter.h"
#include "../Uno_C_API.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// --------------------
// Internal problem struct
// --------------------
struct IpoptProblemInfo {
    ipindex n;
    ipindex m;
    ipindex nele_jac;
    ipindex nele_hess;
    ipindex index_style;

    ipnumber* x_L;
    ipnumber* x_U;
    ipnumber* g_L;
    ipnumber* g_U;

    ipindex* jac_iRow;
    ipindex* jac_jCol;

    ipindex* hess_iRow;
    ipindex* hess_jCol;

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
static uno_int objective_wrapper(uno_int number_variables, const double* x, double* objective_value, void* user_data) {
    IpoptProblem prob = (IpoptProblem)user_data;
    return prob->eval_f((ipindex)number_variables, (const ipnumber*)x, true, objective_value, prob->user_data) ? 0 : 1;
}

static uno_int gradient_wrapper(uno_int number_variables, const double* x, double* gradient, void* user_data) {
    IpoptProblem prob = (IpoptProblem)user_data;
    return prob->eval_grad_f((ipindex)number_variables, (const ipnumber*)x, true, gradient, prob->user_data) ? 0 : 1;
}

static uno_int constraints_wrapper(uno_int number_variables, uno_int number_constraints, const double* x, double* constraint_values, void* user_data) {
    IpoptProblem prob = (IpoptProblem)user_data;
    return prob->eval_g((ipindex)number_variables, (const ipnumber*)x, true,
                        (ipindex)number_constraints, (ipnumber*)constraint_values, prob->user_data) ? 0 : 1;
}

static uno_int jacobian_wrapper(uno_int number_variables, uno_int number_jacobian_nonzeros, const double* x, double* jacobian_values, void* user_data) {
    IpoptProblem prob = (IpoptProblem)user_data;
    return prob->eval_jac_g((ipindex)number_variables, (const ipnumber*)x, true,
                            prob->m, (ipindex)number_jacobian_nonzeros,
                            NULL, NULL, jacobian_values, prob->user_data) ? 0 : 1;
}

static uno_int hessian_wrapper(uno_int number_variables, uno_int number_constraints,
                               uno_int number_hessian_nonzeros,
                               const double* x,
                               double objective_multiplier,
                               const double* multipliers,
                               double* hessian_values,
                               void* user_data)
{
    IpoptProblem prob = (IpoptProblem)user_data;
    return prob->eval_h((ipindex)number_variables, (const ipnumber*)x, true,
        (ipnumber)objective_multiplier, (ipindex)number_constraints, (const ipnumber*)multipliers,
        true, (ipindex)number_hessian_nonzeros, NULL, NULL, hessian_values, prob->user_data) ? 0 : 1;
}

// --------------------
// Free helper
// --------------------
static void free_prob(IpoptProblem p)
{
    if (!p) return;

    free(p->x_L); free(p->x_U);
    free(p->g_L); free(p->g_U);
    free(p->jac_iRow); free(p->jac_jCol);
    free(p->hess_iRow); free(p->hess_jCol);

    if (p->uno_model)  uno_destroy_model(p->uno_model);
    if (p->uno_solver) uno_destroy_solver(p->uno_solver);

    free(p);
}

static uno_int map_index_style(ipindex index_style)
{
    switch(index_style) {
        case 0: return UNO_ZERO_BASED_INDEXING;
        case 1: return UNO_ONE_BASED_INDEXING;
        default: return UNO_ZERO_BASED_INDEXING; // fallback safe
    }
}

// --------------------
// Create
// --------------------
IpoptProblem CreateIpoptProblem(
    ipindex n, const ipnumber* x_L, const ipnumber* x_U,
    ipindex m, const ipnumber* g_L, const ipnumber* g_U,
    ipindex nele_jac, ipindex nele_hess, ipindex index_style,
    Eval_F_CB eval_f, Eval_G_CB eval_g, Eval_Grad_F_CB eval_grad_f,
    Eval_Jac_G_CB eval_jac_g, Eval_H_CB eval_h)
{
    IpoptProblem p = (IpoptProblem)calloc(1, sizeof(struct IpoptProblemInfo));
    if (!p) return NULL;

    p->n = n; p->m = m;
    p->nele_jac = nele_jac;
    p->nele_hess = nele_hess;
    p->index_style = index_style;

    p->eval_f = eval_f;
    p->eval_g = eval_g;
    p->eval_grad_f = eval_grad_f;
    p->eval_jac_g = eval_jac_g;
    p->eval_h = eval_h;

    // copy bounds
    p->x_L = (ipnumber*)malloc(n*sizeof(ipnumber));
    p->x_U = (ipnumber*)malloc(n*sizeof(ipnumber));
    p->g_L = (ipnumber*)malloc(m*sizeof(ipnumber));
    p->g_U = (ipnumber*)malloc(m*sizeof(ipnumber));

    memcpy(p->x_L, x_L, n*sizeof(ipnumber));
    memcpy(p->x_U, x_U, n*sizeof(ipnumber));
    memcpy(p->g_L, g_L, m*sizeof(ipnumber));
    memcpy(p->g_U, g_U, m*sizeof(ipnumber));

    // --------------------
    // Jacobian structure
    // --------------------
    if (nele_jac > 0 && eval_jac_g) {
        p->jac_iRow = (ipindex*)malloc(nele_jac*sizeof(ipindex));
        p->jac_jCol = (ipindex*)malloc(nele_jac*sizeof(ipindex));

        eval_jac_g(n, NULL, true, m, nele_jac,
                   p->jac_iRow, p->jac_jCol,
                   NULL, NULL);

        if (index_style == 1) {
            for (int i=0;i<nele_jac;i++){
                p->jac_iRow[i]--;
                p->jac_jCol[i]--;
            }
        }
    }

    // --------------------
    // Hessian structure
    // --------------------
    if (nele_hess > 0 && eval_h) {
        p->hess_iRow = (ipindex*)malloc(nele_hess*sizeof(ipindex));
        p->hess_jCol = (ipindex*)malloc(nele_hess*sizeof(ipindex));

        eval_h(n, NULL, true,
               1.0, m, NULL, true,
               nele_hess,
               p->hess_iRow, p->hess_jCol,
               NULL, NULL);

        if (index_style == 1) {
            for (int i=0;i<nele_hess;i++){
                p->hess_iRow[i]--;
                p->hess_jCol[i]--;
            }
        }
    }

    // --------------------
    // Create Uno model
    // --------------------
    p->uno_model = uno_create_model(
        UNO_PROBLEM_NONLINEAR,
        n, p->x_L, p->x_U,
        map_index_style(index_style)
    );

    uno_set_objective(p->uno_model,
                      UNO_MINIMIZE,
                      objective_wrapper,
                      gradient_wrapper);

    if (m > 0) {
        uno_set_constraints(
            p->uno_model,
            m,
            constraints_wrapper,
            p->g_L,
            p->g_U,
            nele_jac,
            p->jac_iRow,
            p->jac_jCol,
            jacobian_wrapper
        );
    }

    if (p->eval_h != NULL) {
        uno_set_lagrangian_hessian(
            p->uno_model,
            nele_hess,
            UNO_LOWER_TRIANGLE,
            p->hess_iRow,
            p->hess_jCol,
            hessian_wrapper
        );
    }

    uno_set_lagrangian_sign_convention(p->uno_model, UNO_MULTIPLIER_POSITIVE);

    p->uno_solver = uno_create_solver();
    uno_set_solver_preset(p->uno_solver, "ipopt");

    return p;
}

// --------------------
// Free
// --------------------
void FreeIpoptProblem(IpoptProblem prob)
{
    free_prob(prob);
}

// --------------------
// Options
// --------------------
bool AddIpoptStrOption(IpoptProblem prob, const char* k, const char* v)
{
    return prob && prob->uno_solver &&
           uno_set_solver_string_option(prob->uno_solver, k, v);
}

bool AddIpoptNumOption(IpoptProblem prob, const char* k, ipnumber v)
{
    return prob && prob->uno_solver &&
           uno_set_solver_double_option(prob->uno_solver, k, v);
}

bool AddIpoptIntOption(IpoptProblem prob, const char* k, ipindex v)
{
    return prob && prob->uno_solver &&
           uno_set_solver_integer_option(prob->uno_solver, k, v);
}

bool SetIntermediateCallback(IpoptProblem prob, Intermediate_CB cb)
{
    if (!prob) return false;
    prob->intermediate_cb = cb;
    return true;
}

// --------------------
// Status
// --------------------
ApplicationReturnStatus MapUnoStatus(uno_int termination_status, uno_int solution_status)
{
    // Termination status
    if (termination_status == UNO_ITERATION_LIMIT) {
        return Maximum_Iterations_Exceeded;
    } else if (termination_status == UNO_TIME_LIMIT) {
        return Maximum_Iterations_Exceeded;
    } else if (termination_status == UNO_EVALUATION_ERROR) {
        return Invalid_Problem_Definition;
    } else if (termination_status == UNO_ALGORITHMIC_ERROR) {
        return Internal_Error;
    } else if (termination_status == UNO_USER_TERMINATION) {
        return User_Requested_Stop;
    }

    // Solution status
    if (solution_status == UNO_FEASIBLE_KKT_POINT || solution_status == UNO_FEASIBLE_FJ_POINT) {
        return Solve_Succeeded;
    } else if (solution_status == UNO_INFEASIBLE_STATIONARY_POINT) {
        return Infeasible_Problem_Detected;
    } else if (solution_status == UNO_FEASIBLE_SMALL_STEP || solution_status == UNO_INFEASIBLE_SMALL_STEP) {
        return Search_Direction_Becomes_Too_Small;
    } else if (solution_status == UNO_UNBOUNDED) {
        return Diverging_Iterates;
    } else if (solution_status == UNO_NOT_OPTIMAL) {
        return Unknown_Error;
    } else {
        return Unknown_Error;
    }
}

// --------------------
// Solve
// --------------------

ApplicationReturnStatus IpoptSolve(
    IpoptProblem prob,
    ipnumber* x,
    ipnumber* g,
    ipnumber* obj,
    ipnumber* mult_g,
    ipnumber* mult_x_L,
    ipnumber* mult_x_U,
    UserDataPtr user_data)
{
    if (!prob) return Unknown_Error;

    prob->user_data = user_data;
    uno_set_user_data(prob->uno_model, prob);

    uno_optimize(prob->uno_solver, prob->uno_model);

    // x
    uno_get_primal_solution(prob->uno_solver, x);

    // g (correct)
    if (g && prob->eval_g)
        prob->eval_g(prob->n, x, true, prob->m, g, user_data);

    if (obj)
        *obj = uno_get_solution_objective(prob->uno_solver);

    if (mult_g)
        uno_get_constraint_dual_solution(prob->uno_solver, mult_g);

    if (mult_x_L)
        uno_get_lower_bound_dual_solution(prob->uno_solver, mult_x_L);

    if (mult_x_U)
        uno_get_upper_bound_dual_solution(prob->uno_solver, mult_x_U);

    uno_int term_status = uno_get_optimization_status(prob->uno_solver);
    uno_int sol_status  = uno_get_solution_status(prob->uno_solver);

    return MapUnoStatus(term_status, sol_status);
}

// --------------------
// Version
// --------------------
void GetIpoptVersion(int* major, int* minor, int* release)
{
    if (!major || !minor || !release) return;
    uno_get_version((uno_int*)major,
                    (uno_int*)minor,
                    (uno_int*)release);
}
