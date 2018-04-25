#  NLP written by GAMS Convert at 06/20/02 12:14:12
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
#        22      14       8       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1;
var x3 >= 0, <= 8;
var x4 >= 0;
var x5 >= 0;
var x6 >= 0;
var x7 >= 0;
var x8 >= 0;
var x9 >= 0;

minimize obj: (x3 - 3)*(x3 - 3) + (x1 - 2)*(x1 - 2);

subject to

e2:    x1 - 2*x3 + x4 = 1;

e3:  - 2*x1 + x3 + x5 = 2;

e4:    2*x1 + x3 + x6 = 14;

e5: x4*x7 = 0;

e6: x5*x8 = 0;

e7: x6*x9 = 0;

e8:    2*x1 + x7 - 2*x8 + 2*x9 = 10;
