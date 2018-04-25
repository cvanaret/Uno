#  NLP written by GAMS Convert at 06/20/02 12:22:50
#  
#  Equation counts
#     Total       E       G       L       N       X
#        26      20       3       3       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        30      30       0       0       0       0       0       0
#  FX     2       2       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        81      60      21       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 1000, >= 1000, <= 1000;
var x2 := 1000;
var x3 := 1000;
var x4 := 1000;
var x5 := 100, >= 100;
var x6 := 100, >= 100;
var x7 := 100, >= 100;
var x8 := 400, >= 100, <= 400;
var x9 := 400, >= 100, <= 400;
var x10 := 400, >= 100, <= 400;
var x11 := 100, >= 100, <= 100;
var x12 := 100;
var x13 := 100;
var x14 := 100;
var x15;
var x16;
var x17;
var x19;
var x20;
var x21;
var x22;
var x23;
var x24;
var x25 := 1, >= 0.01;
var x26 := 1, >= 0.01;
var x27 := 1, >= 0.01;
var x28;
var x29;
var x30;

minimize obj:  - x25*x26^0.944*x27^0.891136;

subject to

e2:  - 0.01*(0.5*x5**0.5 + 0.5*(1004.72366 + (-x8) - x15)**0.5)^2 + x25 = 0;

e3:  - 0.01*(0.5*x6**0.5 + 0.5*(1004.72366 + (-x9) - x16)**0.5)^2 + x26 = 0;

e4:  - 0.01*(0.5*x7**0.5 + 0.5*(1004.72366 + (-x10) - x17)**0.5)^2 + x27 = 0;

e5:  - 0.07*x2 - x8 + x28 = 0;

e6:  - 0.07*x3 - x9 + x29 = 0;

e7:  - 0.07*x4 - x10 + x30 = 0;

e8:    x22 - 0.2*x28 = 0;

e9:    x23 - 0.2*x29 = 0;

e10:    x24 - 0.2*x30 = 0;

e11:    x5 + x19 + x22 - x28 = 0;

e12:    x6 + x20 + x23 - x29 = 0;

e13:    x7 + x21 + x24 - x30 = 0;

e14:    x1 - x2 + x11 - x12 + x19 = 0;

e15:    x2 - x3 + x12 - x13 + x20 = 0;

e16:    x3 - x4 + x13 - x14 + x21 = 0;

e17: x15*(x12 - 0.255905*x5) = 1;

e18: x16*(x13 - 0.255905*x6) = 1;

e19: x17*(x14 - 0.255905*x7) = 1;

e20:    x4 + x14 = 1100;

e21:  - 0.25846405*x5 + x12 >= 0;

e22:  - 0.25846405*x6 + x13 >= 0;

e23:  - 0.25846405*x7 + x14 >= 0;

e24:    x8 + x15 <= 904.251294;

e25:    x9 + x16 <= 904.251294;

e26:    x10 + x17 <= 904.251294;
