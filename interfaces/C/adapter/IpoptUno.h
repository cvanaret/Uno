#ifndef IPOPT_UNO_ADAPTER_H
#define IPOPT_UNO_ADAPTER_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef double ipnumber;   // Ipopt Number
typedef int    ipindex;    // Ipopt Index
typedef void*  UserDataPtr;

struct IpoptProblemInfo;
typedef struct IpoptProblemInfo* IpoptProblem;

typedef bool Bool;

// Callback types
typedef bool (*Eval_F_CB)(ipindex n, const ipnumber* x, bool new_x, ipnumber* obj_value, UserDataPtr user_data);
typedef bool (*Eval_Grad_F_CB)(ipindex n, const ipnumber* x, bool new_x, ipnumber* grad_f, UserDataPtr user_data);
typedef bool (*Eval_G_CB)(ipindex n, const ipnumber* x, bool new_x, ipindex m, ipnumber* g, UserDataPtr user_data);
typedef bool (*Eval_Jac_G_CB)(ipindex n, const ipnumber* x, bool new_x, ipindex m, ipindex nele_jac,
                              ipindex* iRow, ipindex* jCol, ipnumber* values, UserDataPtr user_data);
typedef bool (*Eval_H_CB)(ipindex n, const ipnumber* x, bool new_x, ipnumber obj_factor, ipindex m,
                          const ipnumber* lambda, bool new_lambda, ipindex nele_hess,
                          ipindex* iRow, ipindex* jCol, ipnumber* values, UserDataPtr user_data);

typedef bool (*Intermediate_CB)(ipindex alg_mod, ipindex iter_count, ipnumber obj_value,
                                ipnumber inf_pr, ipnumber inf_du, ipnumber mu,
                                ipnumber d_norm, ipnumber regularization_size,
                                ipnumber alpha_du, ipnumber alpha_pr,
                                ipindex ls_trials, UserDataPtr user_data);

// Create / Free
IpoptProblem CreateIpoptProblem(ipindex n, const ipnumber* x_L, const ipnumber* x_U,
                                ipindex m, const ipnumber* g_L, const ipnumber* g_U,
                                ipindex nele_jac, ipindex nele_hess, ipindex index_style,
                                Eval_F_CB eval_f, Eval_G_CB eval_g, Eval_Grad_F_CB eval_grad_f,
                                Eval_Jac_G_CB eval_jac_g, Eval_H_CB eval_h);

void FreeIpoptProblem(IpoptProblem prob);

// Options
bool AddIpoptStrOption(IpoptProblem prob, const char* keyword, const char* val);
bool AddIpoptNumOption(IpoptProblem prob, const char* keyword, ipnumber val);
bool AddIpoptIntOption(IpoptProblem prob, const char* keyword, ipindex val);

// Intermediate callback
bool SetIntermediateCallback(IpoptProblem prob, Intermediate_CB cb);

// Solve
typedef enum {
    Solve_Succeeded = 0,
    Solved_To_Acceptable_Level = 1,
    Infeasible_Problem_Detected = 2,
    Search_Direction_Becomes_Too_Small = 3,
    Diverging_Iterates = 4,
    User_Requested_Stop = 5,
    Feasible_Point_Found = 6,
    Maximum_Iterations_Exceeded = 7,
    Restoration_Failed = 8,
    Error_In_Step_Computation = 9,
    Invalid_Problem_Definition = 10,
    Internal_Error = 11,
    Unknown_Error = 12
} ApplicationReturnStatus;

ApplicationReturnStatus IpoptSolve(IpoptProblem prob, ipnumber* x, ipnumber* g,
                                   ipnumber* obj_val, ipnumber* mult_g,
                                   ipnumber* mult_x_L, ipnumber* mult_x_U,
                                   UserDataPtr user_data);

// Get current iterate / violations
bool GetIpoptCurrentIterate(IpoptProblem prob, bool scaled,
                            ipindex n, ipnumber* x, ipnumber* z_L, ipnumber* z_U,
                            ipindex m, ipnumber* g, ipnumber* lambda);

bool GetIpoptCurrentViolations(IpoptProblem prob, bool scaled,
                               ipindex n, ipnumber* x_L_violation, ipnumber* x_U_violation,
                               ipnumber* compl_x_L, ipnumber* compl_x_U, ipnumber* grad_lag_x,
                               ipindex m, ipnumber* nlp_constraint_violation, ipnumber* compl_g);

// Version
void GetIpoptVersion(int* major, int* minor, int* release);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // IPOPT_UNO_ADAPTER_H
