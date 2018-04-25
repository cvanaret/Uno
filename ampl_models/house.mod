#  NLP written by GAMS Convert at 06/20/02 12:23:28
#  
#  Equation counts
#     Total       E       G       L       N       X
#         9       5       3       1       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         9       9       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        26      17       9       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 30;
var x2;
var x3;
var x4 := 68, >= 40, <= 68;
var x5;
var x6 := 56, >= 56, <= 100;
var x7 <= 3000;
var x8;

minimize obj:  - x7 - x8;

subject to

e1:  - (x1*x2 + x5*x4) + x7 = 0;

e2:  - x1*x3 + x8 = 0;

e4:  - x2 - x5 + x6 = 0;

e5:    x1 - 0.333333333333333*x4 >= 0;

e6:    x1 - 0.5*x4 <= 0;

e7: x2*(x4 - x1) >= 1500;

e8:  - 0.5*x2 + x3 - x5 = 0;

e9:  - 0.5*x2 + x5 >= 0;
