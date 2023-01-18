#include <vector>

extern "C" {
   void nlpqlp_(const int* number_processors, const int* number_constraints, const int* number_equality_constraints, const int* mmax,
         const int* number_variables, const int* nmax, const int* mnn2, double* x /* to check */,
         double* objective_value, double* constraints, double* objective_gradients, double* active_constraint_gradients,
         double* multipliers, double* variable_lower_bounds, double* variable_upper_bounds, double* hessian_approximation, double* cholesky_diagonal,
         double* tolerance, double* QP_tolerance, double* minimum_step_length, int* max_line_search_iterations, int* max_iterations,
         int* max_nonmonotone_iterations, double* restart_hessian_factor, const int* output_level, int* warm_start_mode,
         int* iout, int* ifail, double* wa, int* lwa, int* kwa, const int* lkwa, int* active, const int* lactiv,
         const bool* LQL, void (*)());

   void ql_();
         //void ql_(M, ME, MMAX, N, NMAX, MNN, C, D, A, B, XL, XU, X, U, EPS, MODE, IOUT, FAIL, IPRINT, WAR, LWAR, IWAR, LIWAR )
   void wrapper_(const int* number_processors, const int* number_constraints, const int* number_equality_constraints, const int* mmax,
               const int* number_variables, const int* nmax, const int* mnn2, double* x, double* objective_value, double* constraints,
               double* objective_gradients, double* active_constraint_gradients, double* multipliers, double* variable_lower_bounds,
               double* variable_upper_bounds, double* hessian_approximation, double* cholesky_diagonal,
               double* tolerance, double* QP_tolerance, double* minimum_step_length, int* max_line_search_iterations, int* max_iterations,
               int* max_nonmonotone_iterations, double* restart_hessian_factor, const int* output_level, int* warm_start_mode, int* iout,
               char* ifile, int* ifail, const bool* LQL,
               void (*)(const int* number_constraints, const int* number_equality_constraints, const int* mmax,
                     const int* number_variables, double* objective_value, double* constraints, double* x, int* active, int* ifail),
               void (*)(const int* number_constraints, const int* number_equality_constraints, const int* mmax,
                     const int* number_variables, double* objective_value, double* constraints, double* objective_gradients,
                     double* active_constraint_gradients, double* x, int* active, double* wa)
               );
}

int main() {
   int number_processors = 1;
   int number_constraints = 2;
   int number_equality_constraints = 0;
   int mmax = 0;
   int number_variables = 2;
   int nmax = 0;
   int mnn2 = 0;
   std::vector<double> x(number_variables);
   double objective_value;
   std::vector<double> constraints(number_constraints);
   std::vector<double> objective_gradients(number_variables);
   std::vector<double> active_constraint_gradients(number_equality_constraints*number_variables);
   std::vector<double> multipliers(number_constraints);
   std::vector<double> variable_lower_bounds(number_variables);
   std::vector<double> variable_upper_bounds(number_variables);
   std::vector<double> hessian_approximation(number_variables);
   std::vector<double> cholesky_diagonal(number_variables);
   double tolerance = 1e-6;
   double QP_tolerance = 1e-6;
   double minimum_step_length = 1e-12;
   int max_line_search_iterations = 20;
   int max_iterations = 10000;
   int max_nonmonotone_iterations = 10;
   double restart_hessian_factor = 1.;
   int output_level = 3;
   int warm_start_mode = 0;
   int iout = 0;
   int ifail = 0;
   std::vector<double> wa(number_variables);
   int lwa = 0;
   std::vector<int> kwa(number_variables);
   int lkwa = 0;
   std::vector<int> active(number_constraints);
   int lactiv = 0;
   bool LQL = false;

   nlpqlp_(&number_processors, &number_constraints, &number_equality_constraints, &mmax, &number_variables, &nmax, &mnn2, x.data(), &objective_value,
         constraints.data(), objective_gradients.data(), active_constraint_gradients.data(), multipliers.data(), variable_lower_bounds.data(),
         variable_upper_bounds.data(), hessian_approximation.data(), cholesky_diagonal.data(), &tolerance, &QP_tolerance, &minimum_step_length,
         &max_line_search_iterations, &max_iterations, &max_nonmonotone_iterations, &restart_hessian_factor, &output_level, &warm_start_mode, &iout,
         &ifail, wa.data(), &lwa, kwa.data(), &lkwa, active.data(), &lactiv, &LQL, ql_);
   return 0;
}
