# Optimization sense
const UNO_MINIMIZE = Cint(1)
const UNO_MAXIMIZE = Cint(-1)

# Lagrange multiplier sign convention
const UNO_MULTIPLIER_POSITIVE = Cdouble(1.0)
const UNO_MULTIPLIER_NEGATIVE = Cdouble(-1.0)

# Problem type: 'L' = Linear, 'Q' = Quadratic, 'N' = Nonlinear
const UNO_PROBLEM_LINEAR    = Cchar('L')
const UNO_PROBLEM_QUADRATIC = Cchar('Q')
const UNO_PROBLEM_NONLINEAR = Cchar('N')

# Base indexing style: 0-based (C) or 1-based (Fortran)
const UNO_ZERO_BASED_INDEXING = Cint(0)
const UNO_ONE_BASED_INDEXING  = Cint(1)

# Triangular part: 'L' = lower, 'U' = upper
const UNO_LOWER_TRIANGLE = Cchar('L')
const UNO_UPPER_TRIANGLE = Cchar('U')

# Optimization status
const UNO_SUCCESS = Cint(0)
const UNO_ITERATION_LIMIT = Cint(1)
const UNO_TIME_LIMIT = Cint(2)
const UNO_EVALUATION_ERROR = Cint(3)
const UNO_ALGORITHMIC_ERROR = Cint(4)

# Iterate status
const UNO_NOT_OPTIMAL = Cint(0)
const UNO_FEASIBLE_KKT_POINT = Cint(1)
const UNO_FEASIBLE_FJ_POINT = Cint(2)
const UNO_INFEASIBLE_STATIONARY_POINT = Cint(3)
const UNO_FEASIBLE_SMALL_STEP = Cint(4)
const UNO_INFEASIBLE_SMALL_STEP = Cint(5)
const UNO_UNBOUNDED = Cint(6)
