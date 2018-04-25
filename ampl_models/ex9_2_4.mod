#  NLP written by GAMS Convert at 06/20/02 12:14:02
#  
#  Equation counts
#     Total       E       G       L       N       X
#         8       8       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         9       9       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        19      13       6       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x2;
var x3 >= 0;
var x4 >= 0;
var x5 >= 0;
var x6 >= 0, <= 200;
var x7 >= 0, <= 200;
var x8 >= 0, <= 200;
var x9 >= 0, <= 200;

minimize obj: (0.5*x4 - 1)*(x4 - 2) + (0.5*x5 - 1)*(x5 - 2);

subject to

e2:  - x3 + x4 + x5 = 0;

e3:  - x4 + x6 = 0;

e4:  - x5 + x7 = 0;

e5: x6*x8 = 0;

e6: x7*x9 = 0;

e7:    x2 + x4 - x8 = 0;

e8:    x2 - x9 = -1;
