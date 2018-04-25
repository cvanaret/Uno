#  MIP written by GAMS Convert at 07/08/02 14:14:40
#  
#  Equation counts
#     Total       E       G       L       N       X       C
#        28      16       0      12       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        30      24       6       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        74      74       0       0
# 
#  Reformualtion has removed 1 variable and 1 equation

# modified to a continuous problem


var x2 >= -200, <= 200;
var x3 >= -200, <= 200;
var x4 >= -200, <= 200;
var x5 >= 0;
var x6 >= 0;
var x7 >= 0;
var x8 >= 0;
var x9 >= 0;
var x10 >= 0;
var x11 >= 0;
var x12 >= 0;
var x13 >= 0, <= 200;
var x14 >= 0, <= 200;
var x15 >= 0, <= 200;
var x16 >= 0, <= 200;
var x17 >= 0, <= 200;
var x18 >= 0, <= 200;
var x19 >= 0, <= 200;
var x20 >= 0, <= 200;
var x21 >= 0, <= 200;
var x22 >= 0, <= 200;
var x23 >= 0, <= 200;
var x24 >= 0, <= 200;
var b25 >= 0, <= 1;
var b26 >= 0, <= 1;
var b27 >= 0, <= 1;
var b28 >= 0, <= 1;
var b29 >= 0, <= 1;
var b30 >= 0, <= 1;

minimize obj:  - 8*x5 - 4*x6 + 4*x7 - 40*x8 - 4*x9;

subject to

e2:  - x7 + x8 + x9 + x10 = 1;

e3:    2*x5 - x7 + 2*x8 - 0.5*x9 + x11 = 1;

e4:    2*x6 + 2*x7 - x8 - 0.5*x9 + x12 = 1;

e5:  - x7 + x13 = 0;

e6:  - x8 + x14 = 0;

e7:  - x9 + x15 = 0;

e8:  - x10 + x16 = 0;

e9:  - x11 + x17 = 0;

e10:  - x12 + x18 = 0;

e11:  - x2 - x3 + 2*x4 - x19 = -1;

e12:    x2 + 2*x3 - x4 - x20 = -1;

e13:    x2 - 0.5*x3 - 0.5*x4 - x21 = -2;

e14:    x2 - x22 = 0;

e15:    x3 - x23 = 0;

e16:    x4 - x24 = 0;

e17:    x19 - 400*b25 <= 0;

e18:    x20 - 400*b26 <= 0;

e19:    x21 - 400*b27 <= 0;

e20:    x22 - 400*b28 <= 0;

e21:    x23 - 400*b29 <= 0;

e22:    x24 - 400*b30 <= 0;

e23:    x13 + 400*b25 <= 400;

e24:    x14 + 400*b26 <= 400;

e25:    x15 + 400*b27 <= 400;

e26:    x16 + 400*b28 <= 400;

e27:    x17 + 400*b29 <= 400;

e28:    x18 + 400*b30 <= 400;
