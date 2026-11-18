# Optimization sense
const UNO_MINIMIZE = Cint(1)
const UNO_MAXIMIZE = Cint(-1)

# Lagrange multiplier sign convention
const UNO_MULTIPLIER_POSITIVE = Cint(1)
const UNO_MULTIPLIER_NEGATIVE = Cint(-1)

# Problem type: "LP" = linear problem, "QP" = quadratic problem, "NLP" = nonlinear problem
const UNO_PROBLEM_LINEAR    = "LP"
const UNO_PROBLEM_QUADRATIC = "QP"
const UNO_PROBLEM_NONLINEAR = "NLP"

# Base indexing style: 0-based (C) or 1-based (Fortran)
const UNO_ZERO_BASED_INDEXING = Cint(0)
const UNO_ONE_BASED_INDEXING  = Cint(1)

# Triangular part: 'L' = lower, 'U' = upper
const UNO_LOWER_TRIANGLE = Cchar('L')
const UNO_UPPER_TRIANGLE = Cchar('U')

const UNO_MAX_STRING_LENGTH = Cint(4095)

# Option type: 0 = integer, 1 = double, 2 = boolean, 3 = string, -1 = not found
const UNO_OPTION_TYPE_INTEGER   = Cint(0)
const UNO_OPTION_TYPE_DOUBLE    = Cint(1)
const UNO_OPTION_TYPE_BOOL      = Cint(2)
const UNO_OPTION_TYPE_STRING    = Cint(3)
const UNO_OPTION_TYPE_NOT_FOUND = Cint(-1)

# Status codes returned by the functions that deal with options
const UNO_OPTIONS_SUCCESS          = Cint(0)
const UNO_OPTIONS_INVALID_ARGUMENT = Cint(1)
const UNO_OPTIONS_NOT_FOUND        = Cint(2)
const UNO_OPTIONS_WRONG_TYPE       = Cint(3)
const UNO_OPTIONS_INVALID_VALUE    = Cint(4)
const UNO_OPTIONS_INVALID_FILE     = Cint(5)
const UNO_OPTIONS_ERROR            = Cint(6)
const UNO_OPTIONS_NOT_SET          = Cint(7)

# Status codes returned by the functions that get the result of an optimization
const UNO_RESULT_SUCCESS          = Cint(0)
const UNO_RESULT_INVALID_ARGUMENT = Cint(1)
const UNO_RESULT_NOT_AVAILABLE    = Cint(2)
const UNO_RESULT_INVALID_INDEX    = Cint(3)
const UNO_RESULT_ERROR            = Cint(4)

# Status codes returned by uno_optimize
const UNO_OPTIMIZE_SUCCESS          = Cint(0)
const UNO_OPTIMIZE_INVALID_ARGUMENT = Cint(1)
const UNO_OPTIMIZE_INVALID_MODEL    = Cint(2)
const UNO_OPTIMIZE_INVALID_OPTIONS  = Cint(3)
const UNO_OPTIMIZE_ERROR            = Cint(4)

# Optimization status
const UNO_SUCCESS           = Cint(0)
const UNO_ITERATION_LIMIT   = Cint(1)
const UNO_TIME_LIMIT        = Cint(2)
const UNO_EVALUATION_ERROR  = Cint(3)
const UNO_ALGORITHMIC_ERROR = Cint(4)
const UNO_USER_TERMINATION  = Cint(5)
const UNO_USER_INTERRUPTION = Cint(6)

# Iterate status
const UNO_NOT_OPTIMAL                 = Cint(0)
const UNO_FEASIBLE_KKT_POINT          = Cint(1)
const UNO_FEASIBLE_FJ_POINT           = Cint(2)
const UNO_INFEASIBLE_STATIONARY_POINT = Cint(3)
const UNO_FEASIBLE_SMALL_STEP         = Cint(4)
const UNO_INFEASIBLE_SMALL_STEP       = Cint(5)
const UNO_DIVERGING_ITERATE           = Cint(6)
const UNO_UNBOUNDED_OBJECTIVE         = Cint(7)
