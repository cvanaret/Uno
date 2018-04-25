#  NLP written by GAMS Convert at 06/20/02 11:39:30
#  
#  Equation counts
#     Total       E       G       L       N       X
#         2       1       0       1       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         6       6       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        11       6       5       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 1, >= 0, <= 1;
var x2 := 1, >= 0, <= 1;
var x3 >= 0, <= 1;
var x4 := 1, >= 0, <= 1;
var x5 >= 0, <= 1;

minimize obj: 42*x1 - 0.5*(100*x1*x1 + 100*x2*x2 + 100*x3*x3 + 100*x4*x4 + 100*
              x5*x5) + 44*x2 + 45*x3 + 47*x4 + 47.5*x5;

subject to

e2:    20*x1 + 12*x2 + 11*x3 + 7*x4 + 4*x5 <= 40;
