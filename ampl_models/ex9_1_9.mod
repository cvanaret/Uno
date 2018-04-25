#  MIP written by GAMS Convert at 07/08/02 14:14:41
#  
#  Equation counts
#     Total       E       G       L       N       X       C
#        17       7       0      10       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        18      13       5       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        42      42       0       0
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
var b14 >= 0, <= 1;
var b15 >= 0, <= 1;
var b16 >= 0, <= 1;
var b17 >= 0, <= 1;
var b18 >= 0, <= 1;

minimize obj:    x2 + x3;

subject to

e2:  - x2 - 0.5*x3 + x4 = -2;

e3:  - 0.25*x2 + x3 + x5 = 2;

e4:    x2 + 0.5*x3 + x6 = 8;

e5:    x2 - 2*x3 + x7 = 2;

e6:  - x3 + x8 = 0;

e7:    x9 - 10*b14 <= 0;

e8:    x10 - 10*b15 <= 0;

e9:    x11 - 10*b16 <= 0;

e10:    x12 - 10*b17 <= 0;

e11:    x13 - 10*b18 <= 0;

e12:    x4 + 10*b14 <= 10;

e13:    x5 + 10*b15 <= 10;

e14:    x6 + 10*b16 <= 10;

e15:    x7 + 10*b17 <= 10;

e16:    x8 + 10*b18 <= 10;

e17:  - 0.5*x9 + x10 + 0.5*x11 - 2*x12 - x13 = 1;
