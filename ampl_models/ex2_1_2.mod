#  NLP written by GAMS Convert at 06/20/02 11:40:03
#  
#  Equation counts
#     Total       E       G       L       N       X
#         3       1       0       2       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         7       7       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        15      10       5       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 >= 0, <= 1;
var x2 := 1, >= 0, <= 1;
var x3 >= 0, <= 1;
var x4 := 1, >= 0, <= 1;
var x5 := 1, >= 0, <= 1;
var x6 := 20, >= 0;

minimize obj: (-0.5*(x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5)) - 10.5*x1 - 7.5*x2
               - 3.5*x3 - 2.5*x4 - 1.5*x5 - 10*x6;

subject to

e2:    6*x1 + 3*x2 + 3*x3 + 2*x4 + x5 <= 6.5;

e3:    10*x1 + 10*x3 + x6 <= 20;
