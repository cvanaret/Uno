#include "AMPLModel.hpp"
#include "Logger.hpp"
#include "Utils.hpp"

ASL_pfgh* generate_asl(std::string file_name) {
    SufDecl suffixes[] = {
        /* suffix "uncertain" for variables or constraints */
        {const_cast<char*> (UNCERTAIN_SUFFIX), 0, ASL_Sufkind_var},
        {const_cast<char*> (UNCERTAINTY_SET_SUFFIX), 0, ASL_Sufkind_con}
    };

    ASL_pfgh* asl = (ASL_pfgh*) ASL_alloc(ASL_read_pfgh);
    FILE* nl = jac0dim_ASL((ASL*) asl, const_cast<char*> (file_name.data()), (fint) file_name.size());
    /* indices start at 0 */
    asl->i.Fortran_ = 0;

    // TODO remove macros
    fint n_discrete = nlogv + niv + nlvbi + nlvci + nlvoi;
    if (0 < n_discrete) {
        WARNING << "Ignoring integrality of " << n_discrete << " variables.\n";
        asl->i.need_nl_ = 0;
    }

    /* preallocate initial primal and dual solutions */
    asl->i.X0_ = (double*) M1zapalloc(sizeof (double)*n_var);
    asl->i.pi0_ = (double*) M1zapalloc(sizeof (double)*n_con);

    /* read the file_name.nl file */
    suf_declare_ASL((ASL*) asl, suffixes, sizeof (suffixes) / sizeof (SufDecl));
    pfgh_read_ASL((ASL*) asl, nl, ASL_findgroups);

    return asl;
}

// generate the ASL object and call the private constructor
AMPLModel::AMPLModel(std::string file_name, int fortran_indexing): AMPLModel(file_name, generate_asl(file_name), fortran_indexing) {
}

AMPLModel::AMPLModel(std::string file_name, ASL_pfgh* asl, int fortran_indexing): Problem(file_name, asl->i.n_var_, asl->i.n_con_), variable_uncertain(asl->i.n_var_), constraint_is_uncertainty_set(asl->i.n_con_), asl_(asl), fortran_indexing(fortran_indexing) {
    /* TODO: avoid using implicit AMPL macros */
    keyword keywords[1];
    int size_keywords = sizeof (keywords) / sizeof (keyword);

    Option_Info info = {const_cast<char*> ("Argonot"),
        const_cast<char*> ("Argonot 2018"),
        const_cast<char*> ("argonot_options"),
        keywords, size_keywords};

    //this->asl_ = generate_asl(file_name, &info);
    this->asl_->i.congrd_mode = 0;

    /* dimensions */
    this->objective_sign = (this->asl_->i.objtype_[0] == 1) ? -1. : 1.;

    /* variables */
    this->generate_variables();

    /* objective function */
    this->initialize_objective();

    /* constraints */
    this->generate_constraints();
    set_function_types(file_name, &info);

    /* Lagrangian Hessian */
    this->initialize_lagrangian_hessian();

    /* evaluations counters */
    this->number_eval_objective = 0;
    this->number_eval_constraints = 0;
    this->number_eval_jacobian = 0;
    this->number_eval_hessian = 0;
}

AMPLModel::~AMPLModel() {
    ASL_free((ASL**) & this->asl_);
}

bool is_discrete(ASL_pfgh* asl, int index) {
    return ((asl->i.nlvb_ - asl->i.nlvbi_ <= index && index < asl->i.nlvb_) ||
            (asl->i.nlvc_ - asl->i.nlvci_ <= index && index < asl->i.nlvc_) ||
            (asl->i.nlvo_ - asl->i.nlvoi_ <= index && index < asl->i.nlvo_) ||
            (asl->i.n_var_ - asl->i.niv_ - asl->i.nbv_ <= index && index < asl->i.n_var_));
}

void AMPLModel::generate_variables() {
    SufDesc* uncertain_suffixes = suf_get_ASL((ASL*) this->asl_, UNCERTAIN_SUFFIX, ASL_Sufkind_var);
    
    for (int i = 0; i < this->number_variables; i++) {
        this->variable_name[i] = var_name_ASL((ASL*) this->asl_, i);
        this->variable_discrete[i] = is_discrete(this->asl_, i);
        double lb = (this->asl_->i.LUv_ != NULL) ? this->asl_->i.LUv_[2 * i] : -INFINITY;
        double ub = (this->asl_->i.LUv_ != NULL) ? this->asl_->i.LUv_[2 * i + 1] : INFINITY;
        if (lb == ub) {
            WARNING << "Variable x" << i << " has identical bounds\n";
        }
        this->variables_bounds[i] = {lb, ub};
        this->variable_uncertain[i] = (uncertain_suffixes->u.i != NULL && uncertain_suffixes->u.i[i] == 1);
    }
    this->determine_bounds_types(this->variables_bounds, this->variable_status);
    return;
}

// TODO: fix this duplication!
void AMPLModel::create_objective_variables(ograd* ampl_variables) {
    /* create the dependency pattern as an associative table (variable index, coefficient) */
    ograd* ampl_variables_tmp = ampl_variables;
    while (ampl_variables_tmp != NULL) {
        this->objective_variables[ampl_variables_tmp->varno] = ampl_variables_tmp->coef;
        ampl_variables_tmp = ampl_variables_tmp->next;
    }
    return;
}

void AMPLModel::create_constraint_variables(int j, cgrad* ampl_variables) {
    /* create the dependency pattern as an associative table (variable index, coefficient) */
    cgrad* ampl_variables_tmp = ampl_variables;
    while (ampl_variables_tmp != NULL) {
        this->constraint_variables[j][ampl_variables_tmp->varno] = ampl_variables_tmp->coef;
        ampl_variables_tmp = ampl_variables_tmp->next;
    }
    return;
}

double AMPLModel::objective(std::vector<double>& x) {
    this->number_eval_objective++;
    int nerror = 0;
    double result = this->objective_sign * (*((ASL*) this->asl_)->p.Objval)((ASL*) this->asl_, 0, x.data(), &nerror);
    if (0 < nerror) {
        throw IEEE_FunctionError();
    }
    return result;
}

/* dense gradient */
std::vector<double> AMPLModel::objective_dense_gradient(std::vector<double>& x) {
    std::vector<double> gradient(this->number_variables);
    int nerror = 0;
    /* compute the AMPL gradient (always in dense format) */
    (*((ASL*) this->asl_)->p.Objgrd)((ASL*) this->asl_, 0, x.data(), gradient.data(), &nerror);
    if (0 < nerror) {
        throw IEEE_GradientError();
    }

    /* if maximization, take the opposite */
    if (this->objective_sign < 0.) {
        for (int i = 0; i < this->number_variables; i++) {
            gradient[i] = -gradient[i];
        }
    }
    return gradient;
}

/* sparse gradient */
std::map<int, double> AMPLModel::objective_sparse_gradient(std::vector<double>& x) {
    /* compute the AMPL gradient (always in dense format) */
    std::vector<double> dense_gradient(this->number_variables);
    int nerror = 0;
    (*((ASL*) this->asl_)->p.Objgrd)((ASL*) this->asl_, 0, x.data(), dense_gradient.data(), &nerror);
    if (0 < nerror) {
        throw IEEE_GradientError();
    }

    /* partial derivatives in same order as variables in this->asl_->i.Ograd_[0] */
    std::map<int, double> gradient;
    ograd* ampl_variables_tmp = this->asl_->i.Ograd_[0];
    while (ampl_variables_tmp != NULL) {
        double partial_derivative = dense_gradient[ampl_variables_tmp->varno];
        /* if maximization, take the opposite */
        if (this->objective_sign < 0.) {
            partial_derivative = -partial_derivative;
        }
        gradient[ampl_variables_tmp->varno] = partial_derivative;
        ampl_variables_tmp = ampl_variables_tmp->next;
    }
    return gradient;
}

void AMPLModel::initialize_objective() {
    this->objective_name = obj_name_ASL((ASL*) this->asl_, 0);
    this->create_objective_variables(this->asl_->i.Ograd_[0]);
    return;
}

double AMPLModel::evaluate_constraint(int j, std::vector<double>& x) {
    int nerror = 0;
    double result = (*((ASL*) this->asl_)->p.Conival)((ASL*) this->asl_, j, x.data(), &nerror);
    if (0 < nerror) {
        throw IEEE_FunctionError();
    }
    return result;
}

std::vector<double> AMPLModel::evaluate_constraints(std::vector<double>& x) {
    this->number_eval_constraints++;
    std::vector<double> constraints(this->number_constraints);
    //for (int j = 0; j < this->number_constraints; j++) {
    //    constraints[j] = this->evaluate_constraint(j, x);
    //}
    int nerror = 0;
    (*((ASL*) this->asl_)->p.Conval)((ASL*) this->asl_, x.data(), constraints.data(), &nerror);
    if (0 < nerror) {
        throw IEEE_FunctionError();
    }
    return constraints;
}

/* dense gradient */
std::vector<double> AMPLModel::constraint_dense_gradient(int j, std::vector<double>& x) {
    int congrd_mode_backup = this->asl_->i.congrd_mode;
    this->asl_->i.congrd_mode = 0; // dense computation

    /* compute the AMPL gradient */
    std::vector<double> gradient(this->number_variables);
    int nerror = 0;
    (*((ASL*) this->asl_)->p.Congrd)((ASL*) this->asl_, j, x.data(), gradient.data(), &nerror);
    if (0 < nerror) {
        throw IEEE_FunctionError();
    }

    this->asl_->i.congrd_mode = congrd_mode_backup;

    return gradient;
}

/* sparse gradient */
std::map<int, double> AMPLModel::constraint_sparse_gradient(int j, std::vector<double>& x) {
    int number_variables = this->constraint_variables[j].size(); // <= size(x)
    int congrd_mode_backup = this->asl_->i.congrd_mode;
    this->asl_->i.congrd_mode = 1; // sparse computation

    /* compute the AMPL gradient */
    std::vector<double> ampl_gradient(number_variables);
    int nerror = 0;
    (*((ASL*) this->asl_)->p.Congrd)((ASL*) this->asl_, j, x.data(), ampl_gradient.data(), &nerror);
    if (0 < nerror) {
        throw IEEE_GradientError();
    }

    /* partial derivatives in ampl_gradient in same order as variables in this->asl_->i.Cgrad_[j] */
    std::map<int, double> gradient;
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

std::vector<std::map<int, double> > AMPLModel::constraints_sparse_jacobian(std::vector<double>& x) {
    this->number_eval_jacobian++;
    std::vector<std::map<int, double> > constraints_jacobian(this->number_constraints);
    for (int j = 0; j < this->number_constraints; j++) {
        constraints_jacobian[j] = this->constraint_sparse_gradient(j, x);
    }
    return constraints_jacobian;
}

void AMPLModel::generate_constraints() {
    SufDesc* uncertain_suffixes = suf_get_ASL((ASL*) this->asl_, UNCERTAINTY_SET_SUFFIX, ASL_Sufkind_con);
    
    for (int j = 0; j < this->number_constraints; j++) {
        this->constraint_name[j] = con_name_ASL((ASL*) this->asl_, j);
        this->create_constraint_variables(j, this->asl_->i.Cgrad_[j]);
        double lb = (this->asl_->i.LUrhs_ != NULL) ? this->asl_->i.LUrhs_[2 * j] : -INFINITY;
        double ub = (this->asl_->i.LUrhs_ != NULL) ? this->asl_->i.LUrhs_[2 * j + 1] : INFINITY;
        this->constraint_bounds[j] = {lb, ub};
        this->constraint_is_uncertainty_set[j] = (uncertain_suffixes->u.i != NULL && uncertain_suffixes->u.i[j] == 1);
    }
    this->determine_bounds_types(this->constraint_bounds, this->constraint_status);
    this->determine_constraints();
    return;
}

void AMPLModel::set_function_types(std::string file_name, Option_Info* option_info) {
    /* allocate a temporary ASL to read Hessian sparsity pattern */
    ASL_pfgh* asl = (ASL_pfgh*) ASL_alloc(ASL_read_fg);
    // char* stub = getstops(file_name, option_info);
    //if (file_name == NULL) {
    //	usage_ASL(option_info, 1);
    //}

    FILE* nl = jac0dim(const_cast<char*> (file_name.data()), (fint) file_name.size());
    /* specific read function */
    qp_read_ASL((ASL*) asl, nl, ASL_findgroups);

    fint* rowq;
    fint* colqp;
    double* delsqp;

    /* constraints */
    if (asl->i.n_con_ != this->number_constraints) {
        throw std::length_error("AMPLModel.set_function_types: inconsistent number of constraints");
    }
    this->constraint_type.reserve(this->number_constraints);

    int current_linear_constraint = 0;
    for (int j = 0; j < this->number_constraints; j++) {
        fint qp = nqpcheck_ASL((ASL*) asl, -(j + 1), &rowq, &colqp, &delsqp);

        if (0 < qp) {
            this->constraint_type[j] = QUADRATIC;
        }
        else if (qp == 0) {
            this->constraint_type[j] = LINEAR;
            this->linear_constraints[j] = current_linear_constraint;
            current_linear_constraint++;
        }
        else {
            this->constraint_type[j] = NONLINEAR;
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
    ASL_free((ASL**) & asl);

    return;
}

void AMPLModel::initialize_lagrangian_hessian() {
    /* compute the maximum number of nonzero elements, provided that all multipliers are non-zero */
    /* fint (*Sphset) (ASL*, SputInfo**, int nobj, int ow, int y, int uptri); */
    int obj_number = 0;
    int upper_triangular = 1;
    this->hessian_maximum_number_nonzeros = (*((ASL*) this->asl_)->p.Sphset)((ASL*) this->asl_, NULL, obj_number, 1, 1, upper_triangular);

    /* build sparse description */
    this->hessian_column_start.resize(this->number_variables + 1);
    int* ampl_column_start = this->asl_->i.sputinfo_->hcolstarts;
    for (int k = 0; k < this->number_variables + 1; k++) {
        this->hessian_column_start[k] = ampl_column_start[k] + this->fortran_indexing;
    }

    this->hessian_row_number.resize(this->hessian_maximum_number_nonzeros);
    int* ampl_row_number = this->asl_->i.sputinfo_->hrownos;
    for (int k = 0; k < this->hessian_maximum_number_nonzeros; k++) {
        this->hessian_row_number[k] = ampl_row_number[k] + this->fortran_indexing;
    }

    // use Lagrangian scale: in AMPL, the Lagrangian is f + lambda.g, while Argonot uses f - lambda.g
    int nerror;
    lagscale_ASL((ASL*) this->asl_, -1., &nerror);
    return;
}

CSCMatrix AMPLModel::lagrangian_hessian(std::vector<double>& x, double objective_multiplier, std::vector<double>& multipliers) {
    this->number_eval_hessian++;
    /* register the vector of variables */
    (*((ASL*) this->asl_)->p.Xknown)((ASL*) this->asl_, x.data(), 0);

    /* set the multiplier for the objective function */
    int obj_number = (objective_multiplier != 0.) ? 0 : -1;
    double* obj_multiplier = (objective_multiplier != 0.) ? &objective_multiplier : NULL;

    /* compute the Hessian */
    std::vector<double> hessian(this->hessian_maximum_number_nonzeros);
    (*((ASL*) this->asl_)->p.Sphes)((ASL*) this->asl_, 0, hessian.data(), obj_number, obj_multiplier, multipliers.data());

    /* unregister the vector of variables */
    this->asl_->i.x_known = 0;

    return CSCMatrix(hessian, this->hessian_column_start, this->hessian_row_number, this->fortran_indexing);
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
    std::vector<double> multipliers(this->number_constraints);
    for (int j = 0; j < this->number_constraints; j++) {
        multipliers[j] = ampl_multipliers0[j];
    }
    return multipliers;
}
