#  NLP written by GAMS Convert at 06/20/02 12:13:41
#  
#  Equation counts
#     Total       E       G       L       N       X
#        16      15       0       1       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        17      17       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        45      33      12       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := -8;
var x2 := -8;
var x4 := 1, >= 0, <= 50;
var x5 >= 0, <= 50;
var x6 >= 0, <= 200;
var x7 >= 0, <= 200;
var x8 >= 0, <= 200;
var x9 >= 0, <= 200;
var x10 >= 0, <= 200;
var x11 >= 0, <= 200;
var x12 >= 0, <= 200;
var x13 >= 0, <= 200;
var x14 >= 0, <= 200;
var x15 >= 0, <= 200;
var x16 >= 0, <= 200;
var x17 >= 0, <= 200;

minimize obj:  - 3*x1 - 3*x2 + 2*x4 + 2*x5 - 60;

subject to

e2:    x1 - 2*x2 + x4 + x5 <= 40;

e3:    2*x1 - x4 + x6 = -10;

e4:    2*x2 - x5 + x7 = -10;

e5:  - x1 + x8 = 10;

e6:    x1 + x9 = 20;

e7:  - x2 + x10 = 10;

e8:    x2 + x11 = 20;

e9: x6*x12 = 0;

e10: x7*x13 = 0;

e11: x8*x14 = 0;

e12: x9*x15 = 0;

e13: x10*x16 = 0;

e14: x11*x17 = 0;

e15:    2*x1 - 2*x4 + 2*x12 - x14 + x15 = -40;

e16:    2*x2 - 2*x5 + 2*x13 - x16 + x17 = -40;
