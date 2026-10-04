#include "asl_pfgh.h"
#include "asl2_session.h"

struct Asl2Session {
	ASL_pfgh *asl;
	double *multipliers; /* copy (ASL takes non-const pointers) */
	double *point;
};

Asl2Session* asl2_read(const char* file_name) {
	ASL_pfgh *asl;
	FILE *nl;
	Asl2Session *session;
	char *stub;
	size_t length = strlen(file_name);
	asl = (ASL_pfgh*)ASL_alloc(ASL_read_pfgh);
	stub = (char*)malloc(length + 1);
	strcpy(stub, file_name);
	nl = jac0dim(stub, (fint)length);
	free(stub);
	if (!nl) return 0;
	want_xpi0 = 3;
	pfgh_read(nl, ASL_return_read_err | ASL_findgroups);
	session = (Asl2Session*)malloc(sizeof(Asl2Session));
	session->asl = asl;
	session->multipliers = (double*)calloc(n_con > 0 ? n_con : 1, sizeof(double));
	session->point = (double*)calloc(n_var, sizeof(double));
	return session;
}

void asl2_close(Asl2Session* session) {
	ASL *a = (ASL*)session->asl;
	free(session->multipliers);
	free(session->point);
	ASL_free(&a);
	free(session);
}

#define USE_SESSION ASL_pfgh *asl = session->asl
int asl2_number_variables(const Asl2Session* session) { USE_SESSION; return n_var; }
int asl2_number_constraints(const Asl2Session* session) { USE_SESSION; return n_con; }
size_t asl2_number_jacobian_nonzeros(const Asl2Session* session) { USE_SESSION; return (size_t)nzc; }
void asl2_initial_point(const Asl2Session* session, double* x) {
	USE_SESSION;
	int j;
	for (j = 0; j < n_var; ++j) x[j] = X0 ? X0[j] : 0.;
}

double asl2_objective(Asl2Session* session, const double* x) {
	USE_SESSION;
	fint error = 0;
	return objval(0, (real*)x, &error);
}
void asl2_objective_gradient(Asl2Session* session, const double* x, double* gradient) {
	USE_SESSION;
	fint error = 0;
	objgrd(0, (real*)x, gradient, &error);
}
void asl2_constraints(Asl2Session* session, const double* x, double* constraints) {
	USE_SESSION;
	fint error = 0;
	conval((real*)x, constraints, &error);
}
void asl2_jacobian(Asl2Session* session, const double* x, double* values) {
	USE_SESSION;
	fint error = 0;
	jacval((real*)x, values, &error);
}
void asl2_jacobian_structure(const Asl2Session* session, int* rows, int* columns) {
	USE_SESSION;
	int i;
	cgrad *cg;
	for (i = 0; i < n_con; ++i) {
		for (cg = Cgrad[i]; cg; cg = cg->next) {
			rows[cg->goff] = i;
			columns[cg->goff] = (int)cg->varno;
		}
	}
}
size_t asl2_hessian_setup(Asl2Session* session) {
	USE_SESSION;
	return (size_t)sphsetup(-1, 1, n_con > 0, 1);
}
void asl2_hessian_structure(const Asl2Session* session, int* rows, int* columns) {
	USE_SESSION;
	int j;
	fint k;
	for (j = 0; j < n_var; ++j) {
		for (k = sputinfo->hcolstarts[j]; k < sputinfo->hcolstarts[j + 1]; ++k) {
			rows[k] = (int)sputinfo->hrownos[k];
			columns[k] = j;
		}
	}
}
void asl2_hessian(Asl2Session* session, const double* x, double objective_multiplier, const double* multipliers, double* values) {
	USE_SESSION;
	fint error = 0;
	real ow = objective_multiplier;
	/* ASL evaluates the Hessian at the last point given to a function evaluation */
	xknowne((real*)x, &error);
	if (n_con > 0) memcpy(session->multipliers, multipliers, n_con * sizeof(double));
	sphes(values, -1, &ow, n_con > 0 ? session->multipliers : 0);
	xunknown();
}
void asl2_hessian_vector_product(Asl2Session* session, const double* x, double objective_multiplier,
		const double* multipliers, const double* vector, double* result) {
	USE_SESSION;
	fint error = 0;
	real ow = objective_multiplier;
	xknowne((real*)x, &error);
	if (n_con > 0) memcpy(session->multipliers, multipliers, n_con * sizeof(double));
	hvinit(-1, &ow, n_con > 0 ? session->multipliers : 0);
	hvcomp(result, (real*)vector, -1, &ow, n_con > 0 ? session->multipliers : 0);
	xunknown();
}
