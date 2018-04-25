#  NLP written by GAMS Convert at 06/20/02 11:42:25
#  
#  Equation counts
#     Total       E       G       L       N       X
#         2       2       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        11      11       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        21      11      10       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 >= 0;
var x2 >= 0;
var x3 >= 0;
var x4 := 0.25, >= 0;
var x5 := 0.25, >= 0;
var x6 := 0.25, >= 0;
var x7 := 0.25, >= 0;
var x8 >= 0;
var x9 >= 0;
var x10 >= 0;

minimize obj:  - (x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x6 + x6*x7 + x7*x8 + x8*x9
               + x9*x10 + x1*x3 + x2*x4 + x3*x5 + x4*x6 + x5*x7 + x6*x8 + x7*x9
               + x8*x10 + x1*x9 + x1*x10 + x2*x10 + x1*x5 + x4*x7);

subject to

e2:    x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 = 1;
