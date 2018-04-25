#  MIP written by GAMS Convert at 07/08/02 14:14:40
#  
#  Equation counts
#     Total       E       G       L       N       X       C
#        22      10       0      12       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        24      18       6       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        62      62       0       0
# 
#  Reformualtion has removed 1 variable and 1 equation

# modified to a continuous problem


var x2 >= 0;
var x3 >= 0;
var x4 >= 0;
var x5 >= 0;
var x6 >= 0;
var x7 >= 0;
var x8 >= 0;
var x9 >= 0;
var x10 >= 0;
var x11 >= 0;
var x12 >= 0;
var x13 >= 0;
var x14 >= 0;
var x15 >= 0;
var x16 >= 0;
var x17 >= 0;
var x18 >= 0;
var b19 >= 0, <= 1;
var b20 >= 0, <= 1;
var b21 >= 0, <= 1;
var b22 >= 0, <= 1;
var b23 >= 0, <= 1;
var b24 >= 0, <= 1;

minimize obj:  - 8*x2 - 4*x3 + 4*x4 - 40*x5 + 4*x6;

subject to

e2:  - x4 + x5 + x6 + x7 = 1;

e3:    2*x2 - x4 + 2*x5 - 0.5*x6 + x8 = 1;

e4:    2*x3 + 2*x4 - x5 - 0.5*x6 + x9 = 1;

e5:  - x4 + x10 = 0;

e6:  - x5 + x11 = 0;

e7:  - x6 + x12 = 0;

e8:  - x13 - x14 + 2*x15 - x16 = -1;

e9:    x13 + 2*x14 - x15 - x17 = -1;

e10:    x13 - 0.5*x14 - 0.5*x15 - x18 = -2;

e11:    x13 - 10*b19 <= 0;

e12:    x14 - 10*b20 <= 0;

e13:    x15 - 10*b21 <= 0;

e14:    x16 - 10*b22 <= 0;

e15:    x17 - 10*b23 <= 0;

e16:    x18 - 10*b24 <= 0;

e17:    x7 + 10*b19 <= 10;

e18:    x8 + 10*b20 <= 10;

e19:    x9 + 10*b21 <= 10;

e20:    x10 + 10*b22 <= 10;

e21:    x11 + 10*b23 <= 10;

e22:    x12 + 10*b24 <= 10;
