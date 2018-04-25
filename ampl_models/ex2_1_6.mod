#  NLP written by GAMS Convert at 06/20/02 11:41:25
#  
#  Equation counts
#     Total       E       G       L       N       X
#         6       1       0       5       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        11      11       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        57      47      10       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 1, >= 0, <= 1;
var x2 >= 0, <= 1;
var x3 >= 0, <= 1;
var x4 := 1, >= 0, <= 1;
var x5 := 1, >= 0, <= 1;
var x6 := 1, >= 0, <= 1;
var x7 >= 0, <= 1;
var x8 := 1, >= 0, <= 1;
var x9 := 1, >= 0, <= 1;
var x10 := 1, >= 0, <= 1;

minimize obj: 48*x1 - 0.5*(100*x1*x1 + 100*x2*x2 + 100*x3*x3 + 100*x4*x4 + 100*
              x5*x5 + 100*x6*x6 + 100*x7*x7 + 100*x8*x8 + 100*x9*x9 + 100*x10*
              x10) + 42*x2 + 48*x3 + 45*x4 + 44*x5 + 41*x6 + 47*x7 + 42*x8 + 45
              *x9 + 46*x10;

subject to

e2:  - 2*x1 - 6*x2 - x3 - 3*x5 - 3*x6 - 2*x7 - 6*x8 - 2*x9 - 2*x10 <= -4;

e3:    6*x1 - 5*x2 + 8*x3 - 3*x4 + x6 + 3*x7 + 8*x8 + 9*x9 - 3*x10 <= 22;

e4:  - 5*x1 + 6*x2 + 5*x3 + 3*x4 + 8*x5 - 8*x6 + 9*x7 + 2*x8 - 9*x10 <= -6;

e5:    9*x1 + 5*x2 - 9*x4 + x5 - 8*x6 + 3*x7 - 9*x8 - 9*x9 - 3*x10 <= -23;

e6:  - 8*x1 + 7*x2 - 4*x3 - 5*x4 - 9*x5 + x6 - 7*x7 - x8 + 3*x9 - 2*x10 <= -12;
