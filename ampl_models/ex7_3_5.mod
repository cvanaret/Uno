#  NLP written by GAMS Convert at 06/20/02 11:56:50
#  
#  Equation counts
#     Total       E       G       L       N       X
#        16      12       0       4       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        14      14       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        46      21      25       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1;
var x2;
var x3 >= 0, <= 10;
var x4;
var x5;
var x6;
var x7;
var x8;
var x9;
var x10;
var x11;
var x12;
var x13;

minimize obj:    x4;

subject to

e2: x13*x3^8 - x11*x3^6 + x9*x3^4 - x7*x3^2 + x5 = 0;

e3: x12*x3^6 - x10*x3^4 + x8*x3^2 - x6 = 0;

e4:  - x1 - 0.145*x4 <= -0.175;

e5:    x1 - 0.145*x4 <= 0.175;

e6:  - x2 - 0.15*x4 <= -0.2;

e7:    x2 - 0.15*x4 <= 0.2;

e8:  - 4.53*x1^2 + x5 = 0;

e9:  - (5.28*x1^2 + 0.364*x1) + x6 = 0;

e10:  - (5.72*x1^2*x2 + 1.13*x1^2 + 0.425*x1) + x7 = 0;

e11:  - (6.93*x1^2*x2 + 0.0911*x1) + x8 = 0.00422;

e12:  - (1.45*x1^2*x2 + 0.168*x1*x2) + x9 = 0.000338;

e13:  - (1.56*x1^2*x2^2 + 0.00084*x1^2*x2 + 0.0135*x1*x2) + x10 = 1.35E-5;

e14:  - (0.125*x1^2*x2^2 + 1.68e-5*x1^2*x2 + 0.000539*x1*x2) + x11 = 2.7E-7;

e15:  - (0.005*x1^2*x2^2 + 1.08e-5*x1*x2) + x12 = 0;

e16:  - 0.0001*x1^2*x2^2 + x13 = 0;
