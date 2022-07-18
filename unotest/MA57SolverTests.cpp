// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project root for details.

#include <gtest/gtest.h>

//void test() {
/*
A[0][0] = 2; A[0][1] = 3;
A[1][2] = 4; A[1][4] = 6;
A[2][2] = 1; A[2][3] = 5;
A[4][4] = 1;
 */

//int n = 5;
//COOSymmetricMatrix matrix(n, 0);
//matrix.add_term(2., 0, 0);
//matrix.add_term(3., 0, 1);
//matrix.add_term(4., 1, 2);
//matrix.add_term(6., 1, 4);
//matrix.add_term(1., 2, 2);
//matrix.add_term(5., 2, 3);
//matrix.add_term(1., 4, 4);

/* right hand side */
//std::vector<double> rhs = {8., 45., 31., 15., 17.};

//MA57Solver s;
//MA57Factorization data = s.factorize(matrix);
//std::vector<double> solution = s.solve(matrix, rhs, data);

//std::cout << "solution";
//for (unsigned int k = 0; k < solution.size(); k++) {
//    std::cout << " " << solution[k];
//}
//std::cout << '\n';
//return;
//}



//void test_factorization() {
//    int fortran_indexing = 1;
//    int n = 2;
//    COOSymmetricMatrix coo_matrix(n, fortran_indexing);
//    coo_matrix.add_term(6050.0001, 0, 0);
//    coo_matrix.add_term(-2774, 0, 1);
//    coo_matrix.add_term(1e-4, 1, 1);
//
//    MA57Solver solver;
//    MA57Factorization factorization = solver.factorize(coo_matrix);
//    std::cout << "Dimension: " << coo_matrix.dimension << '\n';
//    std::cout << "Singular ? " << factorization.matrix_is_singular() << '\n';
//    std::cout << "Rank ? " << factorization.rank() << '\n';
//    std::cout << "Negative eigenvalues ? " << factorization.number_negative_eigenvalues() << '\n';
//}