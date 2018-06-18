#include "AMPLModel.hpp"
#include "Logger.hpp"

ASL_pfgh* generate_asl(std::string file_name, Option_Info* option_info) {
	SufDecl suffixes[] = {
		/* suffix "uncertain" for variables or constraints */
		{const_cast<char*>(UNCERTAIN_SUFFIX), 0, ASL_Sufkind_var},
		{const_cast<char*>(UNCERTAINTY_SET_SUFFIX), 0, ASL_Sufkind_con}
	};
	
	ASL_pfgh* asl = (ASL_pfgh*) ASL_alloc(ASL_read_pfgh);
	FILE* nl = jac0dim_ASL((ASL*) asl, const_cast<char*>(file_name.data()), (fint) file_name.size());
	/* indicate that row/col indices are as in Fortran (for bqpd.f) */
	asl->i.Fortran_ = 1;

	// TODO remove macros
	fint n_discrete = nlogv + niv + nlvbi + nlvci + nlvoi;
	if (0 < n_discrete) {
		WARNING << "Ignoring integrality of " << n_discrete << " variables.\n";
		asl->i.need_nl_ = 0;
	}
	
	/* preallocate initial primal and dual solutions */
	asl->i.X0_ = (double*) M1zapalloc(sizeof(double)*n_var);
	asl->i.pi0_ = (double*) M1zapalloc(sizeof(double)*n_con);

	/* read the file_name.nl file */
	suf_declare_ASL((ASL*) asl, suffixes, sizeof(suffixes)/sizeof(SufDecl));
	pfgh_read_ASL((ASL*) asl, nl, ASL_findgroups);
	
	return asl;
}

AMPLModel::AMPLModel(std::string file_name): Problem(file_name) {
	/* avoid using implicit AMPL macros */
	keyword keywords[1];
	int size_keywords = sizeof(keywords)/sizeof(keyword);
	
	Option_Info astros_info = {	const_cast<char*>("astros"),
								const_cast<char*>("ASTROS 2017"),
								const_cast<char*>("astros_options"),
								keywords, size_keywords};
	
	this->asl_ = generate_asl(file_name, &astros_info);
	this->asl_->i.congrd_mode = 0;
	
	/* dimensions */
	this->number_variables = this->asl_->i.n_var_;
	this->number_constraints = this->asl_->i.n_con_;
	this->obj_sign = (this->asl_->i.objtype_[0] == 1) ? -1. : 1.;
	
	/* variables */
	this->generate_variables();
	
	/* objective function */
	this->initialize_objective();
	
	/* constraints */
	this->generate_constraints();
	set_function_types(file_name, &astros_info);
	
	/* fix the Jacobian sparsity pattern */
	this->create_jacobian_sparsity();
	
	/* Lagrangian Hessian */
	this->initialize_lagrangian_hessian();
	
	/* evaluations counters */
	this->number_eval_objective = 0;
	this->number_eval_constraints = 0;
	this->number_eval_hessian = 0;
}

AMPLModel::~AMPLModel() {
	ASL_free((ASL**) &this->asl_);
}

bool is_discrete(ASL_pfgh* asl, int index) {
	return ((asl->i.nlvb_ - asl->i.nlvbi_ <= index && index < asl->i.nlvb_) ||
			(asl->i.nlvc_ - asl->i.nlvci_ <= index && index < asl->i.nlvc_) ||
			(asl->i.nlvo_ - asl->i.nlvoi_ <= index && index < asl->i.nlvo_) ||
			(asl->i.n_var_ - asl->i.niv_ - asl->i.nbv_ <= index && index < asl->i.n_var_));
}

void AMPLModel::generate_variables() {
	SufDesc* uncertain_suffixes = suf_get_ASL((ASL*) this->asl_, UNCERTAIN_SUFFIX, ASL_Sufkind_var);

	this->variable_name.reserve(this->number_variables);
	this->variable_discrete.reserve(this->number_variables);
	this->variable_lb.reserve(this->number_variables);
	this->variable_ub.reserve(this->number_variables);
	this->variable_uncertain.reserve(this->number_variables);
	
	for (int i = 0; i < this->number_variables; i++) {
		this->variable_name.push_back(var_name_ASL((ASL*) this->asl_, i));
		this->variable_discrete[i] = is_discrete(this->asl_, i);
		this->variable_lb[i] = (this->asl_->i.LUv_ != NULL) ? this->asl_->i.LUv_[2*i] : -INFINITY;
		this->variable_ub[i] = (this->asl_->i.LUv_ != NULL) ? this->asl_->i.LUv_[2*i+1] : INFINITY;
		this->variable_uncertain[i] = (uncertain_suffixes->u.i != NULL && uncertain_suffixes->u.i[i] == 1);
	}
	return;
}

// TODO: fix this duplication!
std::map<int,double> create_obj_variables(ograd* ampl_variables) {
	/* create the dependency pattern as an associative table (variable index, coefficient) */
	std::map<int,double> variables;
	
	ograd* ampl_variables_tmp = ampl_variables;
	while (ampl_variables_tmp != NULL) {
		variables[ampl_variables_tmp->varno] = ampl_variables_tmp->coef;
		ampl_variables_tmp = ampl_variables_tmp->next;
	}
	return variables;
}

std::map<int,double> create_cstr_variables(cgrad* ampl_variables) {
	/* create the dependency pattern as an associative table (variable index, coefficient) */
	std::map<int,double> variables;
	
	cgrad* ampl_variables_tmp = ampl_variables;
	while (ampl_variables_tmp != NULL) {
		variables[ampl_variables_tmp->varno] = ampl_variables_tmp->coef;
		ampl_variables_tmp = ampl_variables_tmp->next;
	}
	return variables;
}

double AMPLModel::objective(std::vector<double> x) {
	this->number_eval_objective++;
	int nerror = 0;
	double result = this->obj_sign*(*((ASL*) this->asl_)->p.Objval)((ASL*) this->asl_, 0, x.data(), &nerror);
	if (0 < nerror) {
		throw std::invalid_argument("IEEE error in objective function");
	}
	return result;
}

/* dense gradient */
std::vector<double> AMPLModel::objective_dense_gradient(std::vector<double> x) {
	std::vector<double> gradient(x.size());
	int nerror = 0;
	/* compute the AMPL gradient (always in dense format) */
	(*((ASL*) this->asl_)->p.Objgrd)((ASL*) this->asl_, 0, x.data(), gradient.data(), &nerror);
	if (0 < nerror) {
		throw std::invalid_argument("IEEE error in objective dense gradient");
	}

	/* if maximization, take the opposite */
	if (this->obj_sign < 0.) {
		for (unsigned int i = 0; i < x.size(); i++) {
			gradient[i] = -gradient[i];
		}
	}
	return gradient;
}

/* sparse gradient */
std::map<int,double> AMPLModel::objective_sparse_gradient(std::vector<double> x) {
	/* compute the AMPL gradient (always in dense format) */
	std::vector<double> dense_gradient(x.size());
	int nerror = 0;
	(*((ASL*) this->asl_)->p.Objgrd)((ASL*) this->asl_, 0, x.data(), dense_gradient.data(), &nerror);
	if (0 < nerror) {
		throw std::invalid_argument("IEEE error in objective sparse gradient");
	}

	/* partial derivatives in same order as variables in this->asl_->i.Ograd_[0] */
	std::map<int,double> gradient;
	ograd* ampl_variables_tmp = this->asl_->i.Ograd_[0];
	while (ampl_variables_tmp != NULL) {
		double partial_derivative = dense_gradient[ampl_variables_tmp->varno];
		/* if maximization, take the opposite */
		if (this->obj_sign < 0.) {
			partial_derivative = -partial_derivative;
		}
		gradient[ampl_variables_tmp->varno] = partial_derivative;
		ampl_variables_tmp = ampl_variables_tmp->next;
	}
	return gradient;
}

void AMPLModel::initialize_objective() {
	this->objective_name = obj_name_ASL((ASL*) this->asl_, 0);
	this->objective_variables = create_obj_variables(this->asl_->i.Ograd_[0]);
	return;
}

double AMPLModel::evaluate_constraint(int j, std::vector<double> x) {
	int nerror = 0;
	double result = (*((ASL*) this->asl_)->p.Conival)((ASL*) this->asl_, j, x.data(), &nerror);
	if (0 < nerror) {
		throw std::invalid_argument("IEEE error in constraint function " + j);
	}
	
	return result;
}

std::vector<double> AMPLModel::evaluate_constraints(std::vector<double> x) {
	this->number_eval_constraints++;
	std::vector<double> constraints(this->number_constraints);
	for (int j = 0; j < this->number_constraints; j++) {
		constraints[j] = this->evaluate_constraint(j, x);
	}
	return constraints;
}

/* dense gradient */
std::vector<double> AMPLModel::constraint_dense_gradient(int j, std::vector<double> x) {
	int congrd_mode_backup = this->asl_->i.congrd_mode;
	this->asl_->i.congrd_mode = 0; // dense computation
	
	/* compute the AMPL gradient */
	std::vector<double> gradient(x.size());
	int nerror = 0;
	(*((ASL*) this->asl_)->p.Congrd)((ASL*) this->asl_, j, x.data(), gradient.data(), &nerror);
	if (0 < nerror) {
		throw std::invalid_argument("IEEE error in constraint dense gradient " + j);
	}

	this->asl_->i.congrd_mode = congrd_mode_backup;
	
	return gradient;
}

/* sparse gradient */
std::map<int,double> AMPLModel::constraint_sparse_gradient(int j, std::vector<double> x) {
	int number_variables = this->constraint_variables[j].size(); // <= size(x)
	int congrd_mode_backup = this->asl_->i.congrd_mode;
	this->asl_->i.congrd_mode = 1; // sparse computation
	
	/* compute the AMPL gradient */
	std::vector<double> ampl_gradient(number_variables);
	int nerror = 0;
	(*((ASL*) this->asl_)->p.Congrd)((ASL*) this->asl_, j, x.data(), ampl_gradient.data(), &nerror);
	if (0 < nerror) {
		throw std::invalid_argument("IEEE error in constraint sparse gradient " + j);
	}

	/* partial derivatives in ampl_gradient in same order as variables in this->asl_->i.Cgrad_[j] */
	std::map<int,double> gradient;
	cgrad* ampl_variables_tmp = this->asl_->i.Cgrad_[j];
	int cpt = 0;
	while (ampl_variables_tmp != NULL) {
		/* keep the gradient sparse */
		if (ampl_gradient[cpt] != 0.) {
			gradient[ampl_variables_tmp->varno] = ampl_gradient[cpt];
		}
		ampl_variables_tmp = ampl_variables_tmp->next;
		cpt++;
	}
	
	this->asl_->i.congrd_mode = congrd_mode_backup;

	return gradient;
}

/* dense gradients */
std::vector<std::vector<double> > AMPLModel::constraints_jacobian_dense(std::vector<double> x) {
	/* compute the AMPL gradient */
	std::vector<std::vector<double> > jacobian(this->number_constraints);
	
	for (int j = 0; j < this->number_constraints; j++) {
		jacobian[j] = this->constraint_dense_gradient(j, x);
	}
	
	return jacobian;
}

void AMPLModel::create_jacobian_sparsity() {
	int use_fortran = 1;
	
	/* compute Jacobian sparsity. Initial size = n (dense objective gradient), then grows with the constraints */
	this->jacobian_sparsity.resize(1 + this->number_variables);
	
	/* objective gradient */
	for (int i = 0; i < this->number_variables; i++) {
		this->jacobian_sparsity[i+1] = i + use_fortran;
	}
	
	/* constraint gradients: build Jacobian and sparsity pattern in the same order */
	for (int j = 0; j < this->number_constraints; j++) {
		
		for (std::map<int,double>::iterator it = this->constraint_variables[j].begin(); it != this->constraint_variables[j].end(); it++) {
			int variable_index = it->first;
			this->jacobian_sparsity.push_back(variable_index + use_fortran);
		}
	}
	
	/* Jacobian header */
	this->jacobian_sparsity[0] = jacobian_sparsity.size();
	unsigned int total_size = 1;
	this->jacobian_sparsity.push_back(total_size);
	total_size += this->number_variables;
	this->jacobian_sparsity.push_back(total_size);
	
	for (int j = 0; j < this->number_constraints; j++) {
		total_size += this->constraint_variables[j].size();
		this->jacobian_sparsity.push_back(total_size);
	}
	return;
}

void AMPLModel::generate_constraints() {
	SufDesc* uncertain_suffixes = suf_get_ASL((ASL*) this->asl_, UNCERTAINTY_SET_SUFFIX, ASL_Sufkind_con);
	
	this->constraint_name.reserve(this->number_constraints);
	this->constraint_variables.reserve(this->number_constraints);
	this->constraint_lb.reserve(this->number_constraints);
	this->constraint_ub.reserve(this->number_constraints);
	this->constraint_is_uncertainty_set.reserve(this->number_constraints);
	
	for (int j = 0; j < this->number_constraints; j++) {
		constraint_name.push_back(con_name_ASL((ASL*) this->asl_, j));
		constraint_variables.push_back(create_cstr_variables(this->asl_->i.Cgrad_[j]));
		constraint_lb[j] = (this->asl_->i.LUrhs_ != NULL) ? this->asl_->i.LUrhs_[2*j] : -INFINITY;
		constraint_ub[j] = (this->asl_->i.LUrhs_ != NULL) ? this->asl_->i.LUrhs_[2*j+1] : INFINITY;
		constraint_is_uncertainty_set[j] = (uncertain_suffixes->u.i != NULL && uncertain_suffixes->u.i[j] == 1);
	}
	return;
}

void AMPLModel::set_function_types(std::string file_name, Option_Info* option_info) {
	/* allocate a temporary ASL to read Hessian sparsity pattern */
	ASL_pfgh* asl = (ASL_pfgh*) ASL_alloc(ASL_read_fg);
	// char* stub = getstops(file_name, option_info);
	//if (file_name == NULL) {
	//	usage_ASL(option_info, 1);
	//}
	
	FILE* nl = jac0dim(const_cast<char*>(file_name.data()), (fint) file_name.size());
	/* specific read function */
	qp_read_ASL((ASL*) asl, nl, ASL_findgroups);
	
	fint* rowq;
	fint* colqp;
	double* delsqp;
	
	/* constraints */
	if (asl->i.n_con_ != this->number_constraints) {
		throw std::length_error("AMPLModel.set_function_types: inconsistent number of constraints");
	}
	this->constraints_type.reserve(this->number_constraints);
	
	for (int j = 0; j < this->number_constraints; j++) {
		fint qp = nqpcheck_ASL((ASL*) asl, -(j+1), &rowq, &colqp, &delsqp);
		
		if (0 < qp) {
			this->constraints_type[j] = QUADRATIC;
		}
		else if (qp == 0) {
			this->constraints_type[j] = LINEAR;
		}
		else {
			this->constraints_type[j] = NONLINEAR;
		}
	}
	/* objective function */
	fint qp = nqpcheck_ASL((ASL*) asl, 0, &rowq, &colqp, &delsqp);
	if (0 < qp) {
		this->objective_type = QUADRATIC;
	}
	else if (qp == 0) {
		this->objective_type = LINEAR;
	}
	else {
		this->objective_type = NONLINEAR;
	}
	qp_opify_ASL((ASL*) asl);
	
	/* deallocate memory */
	ASL_free((ASL**) &asl);
	
	return;
}

void AMPLModel::initialize_lagrangian_hessian() {
	/* compute the maximum number of nonzero elements, provided that all multipliers are non-zero */
	/* fint (*Sphset) (ASL*, SputInfo**, int nobj, int ow, int y, int uptri); */
	int obj_number = 0;
	int upper_triangular = 1;
	this->hessian_maximum_number_nonzero = (*((ASL*) this->asl_)->p.Sphset)((ASL*) this->asl_, NULL, obj_number, 1, 1, upper_triangular);
	
	/* build sparse description */
	int use_fortran = 1;
	this->hessian_column_start.resize(this->number_variables+1);
	int* ampl_column_start = this->asl_->i.sputinfo_->hcolstarts;
	for (int k = 0; k < this->number_variables+1; k++) {
		this->hessian_column_start[k] = ampl_column_start[k] - use_fortran;
	}
	
	this->hessian_row_number.resize(this->hessian_maximum_number_nonzero);
	int* ampl_row_number = this->asl_->i.sputinfo_->hrownos;
	for (int k = 0; k < this->hessian_maximum_number_nonzero; k++) {
		this->hessian_row_number[k] = ampl_row_number[k] - use_fortran;
	}
	return;
}

CSCMatrix AMPLModel::lagrangian_hessian(std::vector<double> x, double objective_multiplier, std::vector<double> multipliers) {
	this->number_eval_hessian++;
	/* register the vector of variables */
	(*((ASL*) this->asl_)->p.Xknown)((ASL*) this->asl_, x.data(), 0);

	/* set the multiplier for the objective function */
	int obj_number = (objective_multiplier != 0.) ? 0 : -1;
	double* obj_multiplier = (objective_multiplier != 0.) ? &objective_multiplier : NULL;
	
	/* take the opposite of the multiplier. The reason is that in AMPL, the
	* Lagrangian is f + lambda.g, while Argonot uses f - lambda.g */
	std::vector<double> ampl_multipliers(multipliers.size());
	for (unsigned int j = 0; j < multipliers.size(); j++) {
		ampl_multipliers[j] = -multipliers[j];
	}
	
	/* compute the Hessian */
	std::vector<double> hessian(this->hessian_maximum_number_nonzero);
	(*((ASL*) this->asl_)->p.Sphes)((ASL*) this->asl_, 0, hessian.data(), obj_number, obj_multiplier, ampl_multipliers.data());
	
	/* unregister the vector of variables */
	this->asl_->i.x_known = 0;
	
	return CSCMatrix(hessian, this->hessian_column_start, this->hessian_row_number);
}

/* initial primal point */
std::vector<double> AMPLModel::primal_initial_solution() {
	double* ampl_x0 = this->asl_->i.X0_;
	std::vector<double> x(this->number_variables);
	for (int i = 0; i < this->number_variables; i++) {
		x[i] = ampl_x0[i];
	}
	return x;
}

/* initial dual point */
std::vector<double> AMPLModel::dual_initial_solution() {
	double* ampl_multipliers0 = this->asl_->i.pi0_;
	std::vector<double> multipliers(this->number_variables + this->number_constraints);
	/* n first multipliers (related to bound constraints) are zero */
	for (int j = 0; j < this->number_constraints; j++) {
		multipliers[this->number_variables + j] = ampl_multipliers0[j];
	}
	return multipliers;
}
