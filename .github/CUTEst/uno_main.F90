! THIS VERSION: CUTEST 2.7 - 2026-06-02 AT 11:10 GMT.
! Adapted to the status codes of the Uno C API: the option setters (including
! uno_load_solver_option_file and uno_set_solver_preset) return a status code,
! and uno_optimize returns a logical.

#include "cutest_modules.h"
#include "cutest_routines.h"

!  main interface program from CUTEst to the Uno nonlinear optimization package

!  Nick Gould, Alexis Montoison and Charlie Vanaret, May 2026

  PROGRAM UNO_CUTEST
    USE CUTEST_KINDS_precision
!   USE UNO
    INCLUDE 'uno_c.f90'

!  CUTEst input, output and internal i/o units

    INTEGER ( KIND = ip_ ), PARAMETER :: input = 7
    INTEGER ( KIND = ip_ ), PARAMETER :: out = 6
    INTEGER ( KIND = ip_ ), PARAMETER :: io_buffer = 11
    INTEGER ( KIND = ip_ ) :: output_summary = 47
    CHARACTER ( LEN = 10 ) :: summary_filename = 'Uno.res'
    INTEGER ( KIND = ip_ ) :: last_method = 48
    CHARACTER ( LEN = 10 ) :: last_method_filename = 'Uno.method'
    INTEGER ( KIND = ip_ ) :: output_sol = 49
    CHARACTER ( LEN = 10 ) :: solution_filename = 'Uno.sol'

!  define storage types 

    INTEGER ( KIND = ip_ ) status, i, n, m, l_j, l_h, J_ne, H_ne
    INTEGER ( KIND = uno_int ) :: n_c, m_c, nnz_j_c, nnz_h_c
    INTEGER ( KIND = uno_int ) :: major, minor, patch
    INTEGER ( KIND = uno_int ) :: uno_status
    INTEGER ( KIND = uno_int ) :: optimization_status, iterate_status
    INTEGER ( KIND = uno_int) :: number_iterations
    INTEGER ( KIND = uno_int) :: number_objective_evaluations
    INTEGER ( KIND = uno_int) :: number_constraint_evaluations
    INTEGER ( KIND = uno_int) :: number_objective_gradient_evaluations
    INTEGER ( KIND = uno_int) :: number_jacobian_evaluations
    INTEGER ( KIND = uno_int) :: number_hessian_evaluations
    INTEGER ( KIND = uno_int) :: number_subproblems_solved
    INTEGER ( KIND = uno_int ), PARAMETER ::                                   &
      base_indexing = UNO_ONE_BASED_INDEXING
    INTEGER ( KIND = uno_int ), PARAMETER ::                                   &
      optimization_sense = UNO_MINIMIZE
    INTEGER ( KIND = uno_int ), PARAMETER ::                                   &
      lagrangian_sign_convention = UNO_MULTIPLIER_NEGATIVE
    INTEGER ( KIND = uno_int ), ALLOCATABLE, DIMENSION( : ) ::                 &
      J_row, J_col, H_row, H_col
    REAL ( KIND = c_double ), PARAMETER :: zero = 0.0_c_double
    REAL ( KIND = c_double ) :: solution_objective, solution_primal_feasibility
    REAL ( KIND = c_double ) :: solution_stationarity, solution_complementarity
    REAL ( KIND = c_double ) :: cpu_time, z_i
    REAL ( KIND = c_double ), ALLOCATABLE, DIMENSION( : ) ::                   &
      X_0, X_l, X_u, Y_0, C, C_l, C_u, X, Y, Z_l, Z_u
    LOGICAL :: filexst, same_method, new_header
    LOGICAL ( KIND = c_bool ) :: success
    LOGICAL ( KIND = c_bool ) :: write_solution_to_file
    LOGICAL, ALLOCATABLE, DIMENSION( : ) :: EQUATN, LINEAR
    CHARACTER ( LEN = 7 ) :: tool
    CHARACTER ( LEN = 10 ) :: p_name
    CHARACTER ( LEN = 10 ), ALLOCATABLE, DIMENSION( : ) :: x_names, c_names
    CHARACTER ( LEN = 200 ) :: method_name
    CHARACTER ( : ), ALLOCATABLE :: method_description
    CHARACTER ( LEN = * ), PARAMETER :: options_file_name = "uno.opt"
    CHARACTER ( LEN = * ), PARAMETER :: hessian_model = "exact"
    CHARACTER ( LEN = * ), PARAMETER :: problem_type = UNO_PROBLEM_NONLINEAR
    CHARACTER ( LEN = 1 ), PARAMETER ::                                        &
      hessian_triangular_part = UNO_UPPER_TRIANGLE
    CHARACTER ( LEN = 17 ), DIMENSION( 0 : 6 ), PARAMETER :: OPT_STAT =        &
      [ 'Success          ', 'Iteration limit  ', 'Time limit       ',         &
        'Evaluation error ', 'Algorithmic error', 'User termination ',         &
        'User interruption' ]
    CHARACTER ( LEN = 27 ), DIMENSION( 0 : 7 ), PARAMETER :: ITER_STAT =       &
      [ 'Not optimal                ', 'Feasible KKT point         ',          &
        'Feasible FJ point          ', 'Infeasible stationary point',          &
        'Feasible small step        ', 'Infeasible small step      ',          &
        'Diverging iterate          ', 'Unbounded objective        ' ]
    TYPE ( c_ptr ) :: model, solver
    TYPE ( c_funptr ) :: objective, gradient, constraints, jacobian
    TYPE ( c_funptr ) :: lagrangian_hessian, lagrangian_hessian_operator
    TYPE ( c_funptr ) :: jacobian_operator, jacobian_transposed_operator
    PROCEDURE ( uno_objective_callback ), POINTER ::                           &
      objective_callback => NULL( )
    PROCEDURE ( uno_objective_gradient_callback ), POINTER ::                  &
      objective_gradient_callback => NULL( )
    PROCEDURE ( uno_constraints_callback ), POINTER ::                         &
      constraints_callback => NULL( )
    PROCEDURE ( uno_constraints_jacobian_callback ), POINTER ::                &
      constraints_jacobian_callback => NULL( )
    PROCEDURE ( uno_lagrangian_hessian_callback ), POINTER ::                  &
      lagrangian_hessian_callback => NULL( )
    PROCEDURE ( uno_constraints_jacobian_operator_callback ), POINTER ::       &
      constraints_jacobian_operator_callback => NULL( )
    PROCEDURE ( uno_constraints_jacobian_transposed_operator_callback ),       &
      POINTER :: constraints_jacobian_transposed_operator_callback => NULL( )
    PROCEDURE ( uno_lagrangian_hessian_operator_callback ), POINTER ::         &
      lagrangian_hessian_operator_callback => NULL( )

!  only support specific integer and real types

    IF ( ip_ /= uno_int .OR. c_double /= rpc_) THEN
      WRITE( out,                                                              &
        "( ' Only 32-bit integer and double precision currently supported' )" )
      STOP
    END IF

!  open the relevant data input file

    OPEN( input, FILE = 'OUTSDIF.d', FORM = 'FORMATTED', STATUS = 'OLD' )
    REWIND( input )

!  compute problem dimensions

    tool = 'cdimen '
    CALL CUTEST_cdimen_r( status, input, n, m )
    IF ( status /= 0 ) GO TO 10
    n_c = INT( n, KIND = uno_int )
    m_c = INT( m, KIND = uno_int )

!  allocate space 

    ALLOCATE( X_0( n_c ), Y_0( m_c ), X_l( n_c ), X_u( n_c ),                  &
              C_l( m_c ), C_u( m_c ), X( n_c ), Y( m_c ),                      &
              Z_l( n_c ), Z_u( n_c ), EQUATN( m ), LINEAR( m ), STAT = status )
    IF ( status /= 0 ) THEN
      WRITE( out, "( ' Array allocation failed' )" )
      STOP
    END IF

!  compute the starting points, and the lower and upper bounds on the 
!  variables and constraints

    tool = 'csetup '
    CALL CUTEST_csetup_r( status, input, out, io_buffer, n, m, X_0, X_l, X_u,  &
                          Y_0, C_l, C_u, EQUATN, LINEAR, 0_ip_, 1_ip_, 0_ip_ )
    IF ( status /= 0 ) GO TO 10

!  determine the names of the problem, variables and constraints.

    tool = 'pname  '
    CALL CUTEST_pname_r( status, input, p_name )
    IF ( status /= 0 ) GO TO 10

!  compute the numbers of nonzeros in the constraint Jacobian and Hessian

    tool = 'cdimsj '
    CALL CUTEST_cdimscj_r( status, l_j )
    IF ( status /= 0 ) GO TO 10

    tool = 'cdimsh '
    CALL CUTEST_cdimsh_r( status, l_h )
    IF ( status /= 0 ) GO TO 10
    nnz_j_c = INT( l_j, KIND = uno_int )
    nnz_h_c = INT( l_h, KIND = uno_int )

!  allocate space for matrices

    ALLOCATE( J_row( nnz_j_c ), J_col( nnz_j_c ),                              &
              H_row( nnz_h_c ), H_col( nnz_h_c ), STAT = status )
    IF ( status /= 0 ) THEN
      WRITE( out, "( ' Array allocation failed' )" )
      STOP
    END IF

!  compute the sparsity structure of the Jacobian

    tool = 'csjp   '
    CALL CUTEST_csjp_r( status, J_ne, l_j, J_col, J_row )
    IF ( status /= 0 ) GO TO 10

!  compute the sparsity structure of the Hessian

    tool = 'cshp   '
    CALL CUTEST_cshp_r( status, n, H_ne, l_h, H_row, H_col )
    IF ( status /= 0 ) GO TO 10

!  record the version of Uno used

    CALL uno_get_version( major, minor, patch )
    WRITE( out, "( /, ' Calling Uno version: ', I0, '.', I0, '.', I0 )" )      &
      major, minor, patch

!  callbacks for Uno

    objective_callback => objective_cutest
    objective_gradient_callback  => gradient_cutest
    constraints_callback => constraints_cutest
    constraints_jacobian_callback => jacobian_cutest
    constraints_jacobian_operator_callback => jacobian_operator_cutest
    constraints_jacobian_transposed_operator_callback =>                       &
      jacobian_transposed_operator_cutest
    lagrangian_hessian_callback => lagrangian_hessian_cutest
    lagrangian_hessian_operator_callback => lagrangian_hessian_operator_cutest

    objective = C_FUNLOC( objective_callback )
    gradient = C_FUNLOC( objective_gradient_callback )
    constraints = C_FUNLOC( constraints_callback )
    jacobian = C_FUNLOC( constraints_jacobian_callback )
    jacobian_operator = C_FUNLOC( constraints_jacobian_operator_callback )
    jacobian_transposed_operator                                               &
      = C_FUNLOC( constraints_jacobian_transposed_operator_callback )
    lagrangian_hessian = C_FUNLOC( lagrangian_hessian_callback )
    lagrangian_hessian_operator                                                &
      = C_FUNLOC( lagrangian_hessian_operator_callback )

!  model creation

    model = uno_create_model( problem_type, n_c, X_l, X_u, base_indexing )

    success = uno_set_objective(model, optimization_sense, objective, gradient)
    success                                                                    &
      = uno_set_constraints( model, m_c, constraints, C_l, C_u,                &
                             nnz_j_c, J_row, J_col, jacobian )
    success = uno_set_jacobian_operator(model, jacobian_operator )
    success                                                                    &
       = uno_set_jacobian_transposed_operator( model,                          &
                                               jacobian_transposed_operator )
    success = uno_set_initial_primal_iterate( model, X_0 )
    success = uno_set_initial_dual_iterate( model, Y_0 )

!  solver creation

    solver = uno_create_solver( )

!  load option file

    uno_status = uno_load_solver_option_file( solver, options_file_name )
    CALL CHECK_UNO_STATUS( uno_status, UNO_OPTIONS_SUCCESS,                    &
                           'uno_load_solver_option_file' )

!  how to solve with exact Hessian

    success =                                                                  &
      uno_set_lagrangian_hessian( model, nnz_h_c, hessian_triangular_part,     &
                                  H_row, H_col, lagrangian_hessian )
    uno_status =                                                               &
      uno_set_solver_string_option( solver, "hessian_model", hessian_model )
    CALL CHECK_UNO_STATUS( uno_status, UNO_OPTIONS_SUCCESS,                    &
                           'uno_set_solver_string_option' )
    success =                                                                  &
      uno_set_lagrangian_sign_convention( model, lagrangian_sign_convention )

!  call the optimizer

    success = uno_optimize( solver, model )
    IF ( .NOT. success ) ERROR STOP 'uno_optimize failed'

!  recover method description

    WRITE( out, "( /, ' Solved CUTEst problem ', A, ' with Uno version ',     &
   &  I0, '.', I0, '.', I0, ' using a' )" ) TRIM( p_name ), major, minor, patch
    method_description = uno_get_method_description( solver )
    WRITE( 6, "( 1X, A )" ) method_description

!  recover solution statistics

    optimization_status = uno_get_optimization_status( solver )
    WRITE( out, "( /, ' Optimization status:  ', 11X, A )" )                   &
      OPT_STAT( optimization_status )

    iterate_status = uno_get_solution_status( solver )
    WRITE( out, "( ' Iterate status: ', 17X, A )" ) ITER_STAT( iterate_status )

    solution_objective = uno_get_solution_objective( solver )
    WRITE( out, "( ' Solution objective: ', 11X, ES22.14 )" ) solution_objective

    solution_primal_feasibility = uno_get_solution_primal_feasibility( solver )
    WRITE( out, "( ' Primal feasibility at solution:', ES22.14 )" )            &
      solution_primal_feasibility

    solution_stationarity = uno_get_solution_stationarity( solver )
    WRITE( out, "( ' Stationarity at solution:      ', ES22.14 )" )            &
      solution_stationarity

    solution_complementarity = uno_get_solution_complementarity( solver )
    WRITE( out, "( ' Complementarity at solution:   ', ES22.14 )" )            &
      solution_complementarity

    number_iterations = uno_get_number_iterations( solver )
    WRITE( out, "( ' Number iterations ', 15X, I0 )" ) number_iterations

    number_objective_evaluations                                               &
      = uno_get_number_objective_evaluations( solver )
    WRITE( out, "( ' Objective evaluations: ', 10X, I0 )" )                    &
      number_objective_evaluations                                            

    number_constraint_evaluations                                              &
      = uno_get_number_constraint_evaluations( solver )
    WRITE( out, "( ' Constraints evaluations: ', 8X, I0 )" )                   &
      number_constraint_evaluations

    number_objective_gradient_evaluations                                      &
      = uno_get_number_objective_gradient_evaluations( solver )
    WRITE( out, "( ' Objective gradient evaluations: ', 1X, I0 )" )           &
      number_objective_gradient_evaluations

    number_jacobian_evaluations = uno_get_number_jacobian_evaluations( solver )
    WRITE( out, "( ' Jacobian evaluations: ', 11X, I0 )" )                     &
      number_jacobian_evaluations

    number_hessian_evaluations = uno_get_number_hessian_evaluations( solver )
    WRITE( out, "( ' Hessian evaluations: ', 12X, I0 )" )                      &
      number_hessian_evaluations
    number_subproblems_solved                                                  &
      = uno_get_number_subproblem_solved_evaluations( solver )
    WRITE( out, "( ' Number of subproblems solved: ', 3X, I0 )" )              &
      number_subproblems_solved
    cpu_time = uno_get_cpu_time( solver )
    WRITE( out, "( ' CPU time (secs) ', 13X, F9.3 ) ") cpu_time

    write_solution_to_file                                                     &
      = uno_get_solver_bool_option( solver, "write_solution_to_file" )

!  if required, write the solution to a file

    IF ( write_solution_to_file ) THEN
      INQUIRE( FILE = solution_filename, EXIST = filexst )
      IF ( filexst ) THEN
         OPEN( output_sol, FILE = solution_filename, FORM = 'FORMATTED',       &
               STATUS = 'OLD', IOSTAT = status )
      ELSE
         OPEN( output_sol, FILE = solution_filename, FORM = 'FORMATTED',       &
               STATUS = 'NEW', IOSTAT = status )
      END IF
      IF ( status /= 0 ) THEN
        WRITE( out, "( ' IOSTAT = ', I6, ' when opening file ', A,             &
       &  '. Stopping ' )" ) status, solution_filename
        STOP
      END IF
      REWIND( output_sol )

      ALLOCATE( X_names( n ), STAT = status )
      IF ( status /= 0 ) THEN
        WRITE( out, "( ' Array X_names allocation failed' )" )
        STOP
      END IF

!  record the problem and solver details

      WRITE( output_sol, "( ' Problem: ', A )" ) TRIM( p_name )
      WRITE( output_sol, "( ' Uno version: ',  I0, '.', I0, '.', I0 )" )       &
        major, minor, patch
      WRITE( output_sol, "( ' Uno method: ', A )" )  method_description
      WRITE( output_sol, "( ' Objective value: ', ES22.14 )" )                 &
        solution_objective

!  record the primal and dual solutions

      tool = 'varname'
      CALL CUTEST_varnames( status, n, X_names )
      IF ( status /= 0 ) GO TO 10
      CALL uno_get_primal_solution( solver, X )
      CALL uno_get_lower_bound_dual_solution( solver, Z_l )
      CALL uno_get_upper_bound_dual_solution( solver, Z_u )

      WRITE( output_sol, "( /, ' Solution: ', /,                               &
     & '                                       ',                              &
     & '   <---------------- Bounds ----------------> ', /                     &
     & '       # name              value            ',                         &
     & '       Lower                Upper                  Dual' )" )
      DO i = 1, n
        IF ( Z_l( i ) /= zero ) THEN
          z_i = Z_l( i )
        ELSE IF ( Z_u( i ) /= zero ) THEN
          z_i = Z_u( i )
        ELSE
          z_i = zero
        END IF
        WRITE( output_sol, "( I8, 1X, A10, 4ES22.14 )" )                       &
          i, X_names( i ), X( i ), X_l( i ), X_u( i ), z_i
      END DO
      DEALLOCATE( X_names, STAT = status )

!  record the constraints and their Lagrange multipliers

      IF ( m > 0 ) THEN
        ALLOCATE( C_names( m ), C( m ), STAT = status )
        IF ( status /= 0 ) THEN
          WRITE( out, "( ' Array allocation failed' )" )
          STOP
        END IF
        tool = 'conname'
        CALL CUTEST_connames( status, m, C_names )
        IF ( status /= 0 ) GO TO 10
!       tool = 'ccf    '
!       CALL CUTEST_ccf( status, n, m, X, C )
!       IF ( status /= 0 ) GO TO 10
        call uno_get_solution_constraints( solver, C )
        CALL uno_get_constraint_dual_solution( solver, Y )
        WRITE( output_sol, "( /, ' Constraints: ', /,                          &
       & '                                       ',                            &
       & '   <---------------- Bounds ----------------> ', /                   &
       & '       # name              value            ',                       &
       & '       Lower                Upper               Multiplier' )" )
        DO i = 1, m
          WRITE( output_sol, "( I8, 1X, A10, 4ES22.14 )" )                     &
            i, C_names( i ), C( i ), C_l( i ), C_u( i ), Y( i )
        END DO
        DEALLOCATE( C_names, C, STAT = status )
      END IF
      CLOSE( output_sol )
    END IF

!  append a summary of the results to a file if required

    INQUIRE( FILE = last_method_filename, EXIST = filexst )
    IF ( filexst ) THEN
      OPEN( last_method, FILE = last_method_filename,                          &
            FORM = 'FORMATTED', STATUS = 'OLD', POSITION = 'APPEND',           &
            IOSTAT = status )
    ELSE
      OPEN( last_method, FILE = last_method_filename,                          &
            FORM = 'FORMATTED', STATUS = 'NEW', IOSTAT = status )
    END IF

    IF ( status /= 0 ) THEN
      WRITE( out, "( ' IOSTAT = ', I6, ' when opening file ', A,               &
     & '. Stopping ' )" ) status, last_method_filename
      STOP
    END IF

    IF ( filexst ) THEN
      REWIND( last_method )
      READ( last_method, "( A )" ) method_name
      same_method = TRIM( method_name ) == TRIM( method_description )
      IF ( .NOT. same_method ) THEN
        REWIND( last_method )
        WRITE( last_method, "( A )" ) TRIM( method_description )
      END IF
    ELSE
      same_method = .FALSE.
      WRITE( last_method, "( A )" ) TRIM( method_description )
    END IF
    CLOSE( last_method )

    IF ( status /= 0 ) THEN
      WRITE( out, "( ' IOSTAT = ', I6, ' when opening file ', A,               &
     & '. Stopping ' )" ) status, summary_filename
      STOP
    END IF

    INQUIRE( FILE = summary_filename, EXIST = filexst )
    IF ( filexst ) THEN
      OPEN( output_summary, FILE = summary_filename,                           &
            FORM = 'FORMATTED', STATUS = 'OLD', POSITION = 'APPEND',           &
            IOSTAT = status )
      new_header = .NOT. same_method 
    ELSE
      OPEN( output_summary, FILE = summary_filename,                           &
            FORM = 'FORMATTED', STATUS = 'NEW', IOSTAT = status )
      new_header = .TRUE.
    END IF
    IF ( status /= 0 ) THEN
      WRITE( out, "( ' IOSTAT = ', I6, ' when opening file ', A,               &
     & '. Stopping ' )" ) status, summary_filename
      STOP
    END IF

    IF ( new_header ) THEN
      WRITE( output_summary, "( A )" ) TRIM( method_description )
      WRITE( output_summary, "( 'name              n        m        f   ',    &
     &  '     p_infeas d_infeas c_infeas     iter        time  stat' )" )
    END IF  

    SELECT CASE ( optimization_status )
    CASE ( UNO_SUCCESS )
      WRITE( output_summary,                                                   &
        "( A10, 1X, I8, 1X, I8, ES16.8, 3ES9.1, bn, I9, F12.2, I6 )" )         &
        p_name, n, m, solution_objective, solution_primal_feasibility,         &
        solution_stationarity, solution_complementarity, number_iterations,    &
        cpu_time, optimization_status
    CASE DEFAULT
      WRITE( output_summary,                                                   &
        "( A10, 1X, I8, 1X, I8, ES16.8, 3ES9.1, bn, I9, F12.2, I6 )" )         &
        p_name, n, m, solution_objective, solution_primal_feasibility,         &
        solution_stationarity, solution_complementarity, - number_iterations,  &
        - cpu_time, optimization_status
    END SELECT
    CLOSE( output_summary )

!  cleanup

    CALL uno_destroy_solver( solver )
    CALL uno_destroy_model( model )
    DEALLOCATE( X_0, Y_0, X_l, X_u, C_l, C_u, J_row, J_col, H_row, H_col,      &
                X, Y, Z_l, Z_u, EQUATN, LINEAR, stat = status )
    IF ( status /= 0 ) WRITE( out, "( ' Array deallocation failed' )" )
    STOP

!  record tool failure

 10 CONTINUE
    WRITE( out, "( ' call to CUTEST_', A, ' failed' )" ) tool
    STOP

  CONTAINS

!  stop if a call to Uno did not return the success status

    SUBROUTINE CHECK_UNO_STATUS( uno_status, success_status, function_name )
    INTEGER ( KIND = uno_int ), INTENT( IN ) :: uno_status, success_status
    CHARACTER ( LEN = * ), INTENT( IN ) :: function_name
    IF ( uno_status /= success_status ) THEN
      WRITE( out, "( 1X, A, ' failed with status ', I0 )" )                    &
        function_name, uno_status
      ERROR STOP 'Uno call failed'
    END IF
    END SUBROUTINE CHECK_UNO_STATUS

    INCLUDE 'uno_fortran.f90'

!  objective

    FUNCTION objective_cutest( n, X, f, user_data ) RESULT( res ) BIND( C )
    INTEGER ( KIND = uno_int ), VALUE :: n
    REAL ( KIND = c_double) , INTENT( IN ), DIMENSION( * ) :: X
    REAL ( KIND = c_double ), INTENT( OUT ) :: f
    TYPE( c_ptr ), VALUE :: user_data
    INTEGER ( KIND = uno_int ) :: res
    CALL CUTEST_cifn_r( status, n, 0, X, f )
    res = INT( status, KIND = uno_int )
    RETURN
    END FUNCTION objective_cutest

!  gradient

    FUNCTION gradient_cutest( n, X, G, user_data ) RESULT( res ) BIND( C )
    INTEGER ( KIND = uno_int ), VALUE :: n
    REAL ( KIND = c_double ), INTENT( IN ), DIMENSION( * ) :: X
    REAL ( KIND = c_double ), INTENT( OUT ), DIMENSION( * ) :: G
    TYPE( c_ptr ), VALUE :: user_data
    INTEGER ( KIND = uno_int ) :: res
    CALL CUTEST_cigr_r( status, n, 0, X, G )
    res = INT( status, KIND = uno_int )
    RETURN
    END FUNCTION gradient_cutest

!  constraints

    FUNCTION constraints_cutest( n, m, X, C, user_data ) RESULT( res ) BIND( C )
    INTEGER ( KIND = uno_int ), VALUE :: n, m
    REAL ( KIND = c_double ), INTENT( IN ) , DIMENSION( * ):: X
    REAL ( KIND = c_double ), INTENT( OUT ), DIMENSION( * ) :: C
    TYPE( c_ptr ), VALUE :: user_data
    INTEGER ( KIND = uno_int ) :: res
    CALL CUTEST_ccf_r( status, n, m, X, C )
    res = INT( status, KIND = uno_int )
    RETURN
    END FUNCTION constraints_cutest

!  Jacobian

    FUNCTION jacobian_cutest( n, nnz_j, X, J_val,                              &
                              user_data ) RESULT( res ) BIND( C )
    INTEGER ( KIND = uno_int ), VALUE :: n, nnz_j
    REAL ( KIND = c_double ), INTENT( IN ), DIMENSION( * ) :: X
    REAL ( KIND = c_double ), INTENT( OUT ), DIMENSION( * ) :: J_val
    TYPE( c_ptr ), VALUE :: user_data
    INTEGER ( KIND = uno_int ) :: res
    INTEGER ( KIND = ip_ ) :: l_j
    INTEGER ( KIND = ip_ ), DIMENSION( nnz_j ) :: J_row, J_col
    l_j =  nnz_j
    CALL CUTEST_csj_r( status, n, X, nnz_j, l_j, J_val, J_col, J_row )
    res = INT( status, KIND = uno_int )
    RETURN
    END FUNCTION jacobian_cutest

! Jacobian operator

    FUNCTION jacobian_operator_cutest( n, m, X, evaluate_at_x, VECTOR, RESULT, &
                                       user_data ) RESULT( res ) BIND( C )
    INTEGER ( KIND = uno_int ), VALUE :: n, m
    REAL ( KIND = c_double ), INTENT( IN ), DIMENSION( * ) :: X, VECTOR
    LOGICAL ( KIND = c_bool ), VALUE :: evaluate_at_x
    REAL ( KIND = c_double ), INTENT( OUT ), DIMENSION( * ) :: RESULT
    TYPE( c_ptr ), VALUE :: user_data
    INTEGER ( KIND = uno_int ) :: res
    LOGICAL :: got_j
    got_j = .NOT. evaluate_at_x
    CALL CUTEST_cjprod_r( status, n, m, got_j, .FALSE., X, VECTOR, n, RESULT, n)
    res = INT( status, KIND = uno_int )
    RETURN
    END FUNCTION jacobian_operator_cutest

! Jacobian transposed operator

    FUNCTION jacobian_transposed_operator_cutest( n, m, X, evaluate_at_x,      &
                          VECTOR, RESULT, user_data ) RESULT( res ) BIND( C )
    INTEGER ( KIND = uno_int ), VALUE :: n, m
    REAL ( KIND = c_double ), INTENT( IN ), DIMENSION( * ) :: X, VECTOR
   LOGICAL ( KIND = c_bool ), VALUE :: evaluate_at_x
    REAL ( KIND = c_double ), INTENT( OUT ), DIMENSION( * ) :: RESULT
    TYPE( c_ptr ), VALUE :: user_data
    INTEGER ( KIND = uno_int ) :: res
    LOGICAL :: got_j
    got_j = .NOT. evaluate_at_x
    CALL CUTEST_cjprod_r( status, n, m, got_j, .TRUE., X, VECTOR, n, RESULT, n )
    res = INT( status, KIND = uno_int )
    RETURN
    END FUNCTION jacobian_transposed_operator_cutest

!  Lagrangian Hessian

    FUNCTION lagrangian_hessian_cutest( n, m, nnz_h, X, y0, Y, H_val,          &
                                        user_data ) RESULT( res ) BIND( C )
    INTEGER ( KIND = uno_int ), VALUE :: n, m, nnz_h
    REAL ( KIND = c_double ), INTENT( IN ), DIMENSION( * ) :: X, Y
    REAL ( KIND = c_double ), INTENT( OUT ), DIMENSION( * ) :: H_val
    REAL ( KIND = c_double ), VALUE :: y0
    TYPE( c_ptr ), VALUE :: user_data
    INTEGER ( KIND = uno_int ) :: res
    INTEGER ( KIND = ip_ ) :: l_h
    INTEGER ( KIND = ip_ ), DIMENSION( nnz_h ) :: H_row, H_col
    REAL ( KIND = c_double ) :: Y_null( 1 )
    REAL ( KIND = c_double ), DIMENSION( m ) :: YY
    l_h =  nnz_h
    IF ( m > 0 ) THEN
      YY( 1 : m ) = Y( 1 : m )
      CALL CUTEST_cshj_r( status, n, m, X, y0, - YY, nnz_h, l_h,               &
                          H_val, H_row, H_col )
    ELSE
      CALL CUTEST_cshj_r( status, n, m, X, y0, Y_null( 1 ), nnz_h, l_h,        &
                          H_val, H_row, H_col )
    END IF
    res = INT( status, KIND = uno_int )
    RETURN
    END FUNCTION lagrangian_hessian_cutest

!  Lagrangian Hessian operator

    FUNCTION lagrangian_hessian_operator_cutest( n, m, X, evaluate_at_x, y0,   &
                          Y, VECTOR, RESULT, user_data ) RESULT( res ) BIND( C )
    INTEGER ( KIND = uno_int ), VALUE :: n, m
    REAL ( KIND = c_double), INTENT( IN ), DIMENSION( * ) :: X, Y, VECTOR
    LOGICAL ( KIND = c_bool ), VALUE :: evaluate_at_x
    REAL ( KIND = c_double ), VALUE :: y0
    REAL ( KIND = c_double ), INTENT( OUT ), DIMENSION( * ) :: RESULT
    TYPE( c_ptr ), VALUE :: user_data
    INTEGER ( KIND = uno_int ) :: res
    REAL ( KIND = c_double ) :: Y_null( 1 )
    REAL ( KIND = c_double ), DIMENSION( m ) :: YY
    LOGICAL :: got_h
    got_h = .NOT. evaluate_at_x
    IF ( m > 0 ) THEN
      YY( 1 : m ) = Y( 1 : m )
      CALL CUTEST_chjprod_r( status, n, m, got_h, X, y0, - YY,                 &
                             VECTOR, RESULT )
    ELSE
      CALL CUTEST_chjprod_r( status, n, m, got_h, X, y0, Y_null( 1 ),          &
                             VECTOR, RESULT )
    END IF
    res = INT( status, KIND = uno_int )
    RETURN
    END FUNCTION lagrangian_hessian_operator_cutest

  END PROGRAM UNO_CUTEST
