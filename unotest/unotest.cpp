// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#if defined(HAS_MPI) && defined(MUMPS_PARALLEL)
#include "mpi.h"
#endif
#include <gtest/gtest.h>
#include "tools/Logger.hpp"

// https://www.eriksmistad.no/getting-started-with-google-test-on-ubuntu/
int main(int argc, char **argv) {
#if defined(HAS_MPI) && defined(MUMPS_PARALLEL)
   int myid , ierr;
   ierr = MPI_Init(&argc , &argv);
   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#endif

    testing::InitGoogleTest(&argc, argv);
    auto result = RUN_ALL_TESTS();

#if defined(HAS_MPI) && defined(MUMPS_PARALLEL)
   ierr = MPI_Finalize() ;
#endif
   return result;
}

//void test_mask_matrix() {
//    int n = 4;
//    UnoMatrix matrix(n, 0);
//
//    // Column 0
//    matrix.insert(1., 0, 0);
//    // Column 1
//    matrix.insert(2., 0, 1);
//    matrix.insert(3., 1, 1);
//    // Column 2
//    matrix.insert(4., 0, 2);
//    matrix.insert(5., 1, 2);
//    matrix.insert(6., 2, 2);
//    // Column 3
//    matrix.insert(7., 0, 3);
//    matrix.insert(8., 1, 3);
//    matrix.insert(9., 2, 3);
//    matrix.insert(10., 3, 3);
//    std::cout << "Original matrix: " << matrix << '\n';
//
//    std::unordered_map<int, int> mask;
//    mask[0] = 0;
//    mask[2] = 1;
//
//
//    COOSymmetricMatrix coo_matrix = matrix.to_COO(mask);
//    std::cout << "COO reduced matrix:\n" << coo_matrix;
//
//    CSCSymmetricMatrix csc_matrix = matrix.to_CSC(mask);
//    std::cout << "CSC reduced matrix:\n" << csc_matrix;
//}

// C++ problem

//double f(std::vector<double> x) {
//    return x[0];
//}
//
//std::vector<double> f_gradient(std::vector<double> x) {
//    std::vector<double> gradient(1);
//    gradient[0] = 1.;
//    return gradient;
//}

//void test_cpp() {
//    CppProblem problem("my_problem", 1, 0, f, f_gradient);
//    std::vector<double> x = {123.};
//    double f_x = f(x);
//    std::cout << "f(x) = " << f_x << '\n';
//}

// sparse vector

//void dense_to_sparse(const std::vector<double>& input) {
   //// create
   //std::vector<std::pair<int, double>> output;
   //output.reserve(input.size());
   //for (size_t i = 0; i < input.size(); i++) {
      //if (input[i] != 0.) {
         //output.emplace_back(i, input[i]);
      //}
   //}
   //// print
   //for (auto& [i, xi]: output) {
      //std::cout << "Index " << i << ", value " << xi << '\n';
   //}
//}

//void test_sparse_vector() {
   //std::vector<double> input{0, 1, 2, 0, 0, 0, 3};
   //dense_to_sparse(input);
//}


//double f(const std::vector<double>& x) {
    //return 2*x[0];
//}

//void f_gradient(const std::vector<double>& /*x*/, std::vector<double>& gradient) {
    //gradient[0] = 1.;
//}

//void test_cpp() {
    //CppProblem<1, 0, 0> problem("my_problem", f, f_gradient);
    //std::vector<double> x = {123.};
    //double f_x = f(x);
    //std::cout << "f(x) = " << f_x << '\n';
//}

//void test_product() {
//std::vector<double> m(5);
//m[0] = 1.;
//m[1] = 3.;
//m[2] = 4.;
//m[3] = 7.;
//m[4] = 15.;
//std::vector<int> matrix_column_start = {1, 2, 4, 6};
//std::vector<int> matrix_row_number = {1, 1, 2, 1, 3};
//Matrix matrix(m, matrix_column_start, matrix_row_number);
//std::vector<double> vector(3);
//vector[0] = 1.;
//vector[1] = 2.;
//vector[2] = 3.;

//std::vector<double> result = matrix.product(vector);
///* *expected result = (28, 11, 52) */
//}


/*
void test_matrix() {
    int size = 3;
    int nnz = 4;
    UnoMatrix matrix(size, nnz);
    // add terms
    matrix.add_term(12., 0, 0);
    
    matrix.add_term(14., 0, 1);
    matrix.add_term(7., 0, 1);
    
    matrix.add_term(13., 1, 2);
    
    matrix.add_term(19., 2, 2);
    
    print_vector(std::cout, matrix.matrix);
    
    COOSymmetricMatrix coo_matrix = matrix.to_COO();
    std::cout << "COO matrix: \n";
    print_vector(std::cout, coo_matrix.matrix);
    print_vector(std::cout, coo_matrix.row_indices);
    print_vector(std::cout, coo_matrix.column_indices);
}
 */
