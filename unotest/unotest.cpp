// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#if defined(HAS_MPI) && defined(MUMPS_PARALLEL)
#include "mpi.h"
#endif
#include <gtest/gtest.h>

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