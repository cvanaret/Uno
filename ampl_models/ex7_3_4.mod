#  NLP written by GAMS Convert at 06/20/02 11:56:19
#  
#  Equation counts
#     Total       E       G       L       N       X
#        18       8       0      10       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        13      13       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        53      30      23       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1;
var x2;
var x3;
var x4;
var x5;
var x6;
var x7;
var x8;
var x9;
var x10;
var x11 >= 0, <= 10;
var x12;

minimize obj:    x12;

subject to

e2: x10*x11^4 - x8*x11^2 + x6 = 0;

e3: x9*x11^2 - x7 = 0;

e4:  - x1 - x12 <= -10;

e5:    x1 - x12 <= 10;

e6:    x2 - 0.1*x12 <= 1;

e7:  - x2 - 0.1*x12 <= -1;

e8:  - x3 - 0.1*x12 <= -1;

e9:    x3 - 0.1*x12 <= 1;

e10:  - x4 - 0.01*x12 <= -0.2;

e11:    x4 - 0.01*x12 <= 0.2;

e12:  - x5 - 0.005*x12 <= -0.05;

e13:    x5 - 0.005*x12 <= 0.05;

e14:  - 54.387*x3*x2 + x6 = 0;

e15:  - 0.2*(1364.67*x3*x2 - 147.15*x4*x3*x2) + 5.544*x5 + x7 = 0;

e16:  - 3*(-9.81*x3*x2^2 - 9.81*x3*x1*x2 - 4.312*x3^2*x2 + 264.896*x3*x2 + x4*
     x5 - 9.274*x5) + x8 = 0;

e17:  - (7*x4*x3^2*x2 - 64.918*x3^2*x2 + 380.067*x3*x2 + 3*x5*x2 + 3*x5*x1)
      + x9 = 0;

e18:  - x3^2*x2*(7*x1 + 4*x2) + x10 = 0;
