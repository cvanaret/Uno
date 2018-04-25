#  NLP written by GAMS Convert at 06/20/02 11:48:23
#  
#  Equation counts
#     Total       E       G       L       N       X
#         7       5       0       2       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        10      10       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        30      23       7       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 >= 0, <= 100;
var x2 >= 0, <= 200;
var x3 >= 0, <= 500;
var x4 >= 0, <= 500;
var x5 >= 0, <= 500;
var x6 >= 0, <= 500;
var x7 >= 0, <= 500;
var x8 >= 0, <= 500;
var x9 >= 0, <= 500;

minimize obj:  - 9*x1 - 15*x2 + 6*x3 + 13*x4 + 10*x5 + 10*x6;

subject to

e2:  - x3 - x4 + x8 + x9 = 0;

e3:    x1 - x5 - x8 = 0;

e4:    x2 - x6 - x9 = 0;

e5: x7*x8 - 2.5*x1 + 2*x5 <= 0;

e6: x7*x9 - 1.5*x2 + 2*x6 <= 0;

e7: x7*x8 + x7*x9 - 3*x3 - x4 = 0;
