#ifndef CPPASL_ASL2_SESSION_H
#define CPPASL_ASL2_SESSION_H
/* Thin C wrapper around the original ASL2 (solvers2) evaluation API, compiled in its own translation unit so that
 * ASL's macros do not leak into C++ code. */
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct Asl2Session Asl2Session;
/* reads the .nl file with pfgh_read (partially separable, Hessian support); NULL on failure */
Asl2Session* asl2_read(const char* file_name);
void asl2_close(Asl2Session* session);
int asl2_number_variables(const Asl2Session* session);
int asl2_number_constraints(const Asl2Session* session);
size_t asl2_number_jacobian_nonzeros(const Asl2Session* session);
void asl2_initial_point(const Asl2Session* session, double* x);
double asl2_objective(Asl2Session* session, const double* x);
void asl2_objective_gradient(Asl2Session* session, const double* x, double* gradient);
void asl2_constraints(Asl2Session* session, const double* x, double* constraints);
void asl2_jacobian(Asl2Session* session, const double* x, double* values);
/* (row, column) of each Jacobian value, in ASL's (goff) order */
void asl2_jacobian_structure(const Asl2Session* session, int* rows, int* columns);
/* sphsetup for the Lagrangian (objective 0 + all constraints), upper triangle; returns the number of nonzeros */
size_t asl2_hessian_setup(Asl2Session* session);
/* (row, column) with row <= column */
void asl2_hessian_structure(const Asl2Session* session, int* rows, int* columns);
void asl2_hessian(Asl2Session* session, const double* x, double objective_multiplier, const double* multipliers, double* values);
/* Lagrangian Hessian-vector product at x (hvinit + hvcomp) */
void asl2_hessian_vector_product(Asl2Session* session, const double* x, double objective_multiplier,
   const double* multipliers, const double* vector, double* result);
#ifdef __cplusplus
}
#endif
#endif
