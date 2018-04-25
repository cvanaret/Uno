#  NLP written by GAMS Convert at 06/20/02 11:49:28
#  
#  Equation counts
#     Total       E       G       L       N       X
#        17      17       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        23      23       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        64      40      24       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 >= 0, <= 300;
var x2 >= 0, <= 300;
var x3 >= 0, <= 300;
var x4 >= 0, <= 300;
var x5 >= 0, <= 300;
var x6 >= 0, <= 300;
var x7 >= 0, <= 300;
var x8 >= 0, <= 300;
var x9 >= 0, <= 300;
var x10 >= 0, <= 300;
var x11 >= 0, <= 300;
var x12 >= 0, <= 300;
var x13 >= 0, <= 300;
var x14 >= 0, <= 300;
var x15 >= 0, <= 300;
var x16 >= 0, <= 300;
var x17 >= 0, <= 300;
var x18 >= 0, <= 300;
var x19 >= 0, <= 1;
var x20 >= 0, <= 1;
var x21 >= 0, <= 1;
var x22 >= 0, <= 1;

minimize obj:    0.00432*x1 + 0.01517*x2 + 0.01517*x9 + 0.00432*x13 + 0.9979;

subject to

e1:    x1 + x2 + x3 + x4 = 300;

e2:    x5 - x6 - x7 = 0;

e3:    x8 - x9 - x10 - x11 = 0;

e4:    x12 - x13 - x14 - x15 = 0;

e5:    x16 - x17 - x18 = 0;

e6: x13*x21 + 0.333*x1 - x5 = 0;

e7: x13*x22 - x8*x20 + 0.333*x1 = 0;

e8:  - x8*x19 + 0.333*x1 = 0;

e9:  - x12*x21 - 0.333*x2 = 0;

e10: x9*x20 - x12*x22 + 0.333*x2 = 0;

e11: x9*x19 + 0.333*x2 - x16 = 0;

e12: x14*x21 + 0.333*x3 + x6 = 30;

e13: x10*x20 + x14*x22 + 0.333*x3 = 50;

e14: x10*x19 + 0.333*x3 + x17 = 30;

e15:    x19 + x20 = 1;

e16:    x21 + x22 = 1;
