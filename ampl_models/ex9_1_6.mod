#  MIP written by GAMS Convert at 07/08/02 14:14:40
#  
#  Equation counts
#     Total       E       G       L       N       X       C
#        20       8       0      12       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        21      15       6       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        50      50       0       0
# 
#  Reformualtion has removed 1 variable and 1 equation

# modified to a continuous problem


var x2 >= 0;
var x3 >= 0;
var x4 >= 0, <= 200;
var x5 >= 0, <= 200;
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
var b16 >= 0, <= 1;
var b17 >= 0, <= 1;
var b18 >= 0, <= 1;
var b19 >= 0, <= 1;
var b20 >= 0, <= 1;
var b21 >= 0, <= 1;

minimize obj:  - x2 - 3*x3;

subject to

e2:  - x2 - 2*x3 + x4 = -10;

e3:    x2 - 2*x3 + x5 = 6;

e4:    2*x2 - x3 + x6 = 21;

e5:    x2 + 2*x3 + x7 = 38;

e6:  - x2 + 2*x3 + x8 = 18;

e7:  - x3 + x9 = 0;

e8:  - 2*x10 - 2*x11 - x12 + 2*x13 + 2*x14 - x15 = -3;

e9:    x10 - 400*b16 <= 0;

e10:    x11 - 400*b17 <= 0;

e11:    x12 - 400*b18 <= 0;

e12:    x13 - 400*b19 <= 0;

e13:    x14 - 400*b20 <= 0;

e14:    x15 - 400*b21 <= 0;

e15:    x4 + 400*b16 <= 400;

e16:    x5 + 400*b17 <= 400;

e17:    x6 + 400*b18 <= 400;

e18:    x7 + 400*b19 <= 400;

e19:    x8 + 400*b20 <= 400;

e20:    x9 + 400*b21 <= 400;
