#  NLP written by GAMS Convert at 07/06/02 13:18:05
#  
#  Equation counts
#     Total       E       G       L       N       X       C
#         7       2       0       5       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         8       8       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        28      12      16       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 >= 0, <= 1;
var x2 >= 0, <= 1;
var x3 >= 0, <= 1;
var x4 >= 0, <= 100;
var x5 >= 0, <= 200;
var x6 >= 0, <= 100;
var x7 >= 0, <= 200;

minimize obj:  - ((9 + (-6*x1) - 16*x2 - 15*x3)*x4 + (15 + (-6*x1) - 16*x2 - 15
              *x3)*x5) + x6 - 5*x7;

subject to

e2: x3*x4 + x3*x5 <= 50;

e3:    x4 + x6 <= 100;

e4:    x5 + x7 <= 200;

e5: (3*x1 + x2 + x3 - 2.5)*x4 - 0.5*x6 <= 0;

e6: (3*x1 + x2 + x3 - 1.5)*x5 + 0.5*x7 <= 0;

e7:    x1 + x2 + x3 = 1;
