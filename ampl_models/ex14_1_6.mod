#  NLP written by GAMS Convert at 06/20/02 11:34:49
#  
#  Equation counts
#     Total       E       G       L       N       X
#        16       2       0      14       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        10      10       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        62      30      32       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 >= -1, <= 1;
var x2 >= -1, <= 1;
var x3 >= -1, <= 1;
var x4 >= -1, <= 1;
var x5 >= -1, <= 1;
var x6 >= -1, <= 1;
var x7 >= -1, <= 1;
var x8 >= -1, <= 1;
var x9;

minimize obj:    x9;

subject to

e2: 0.004731*x1*x3 - 0.1238*x1 - 0.3578*x2*x3 - 0.001637*x2 - 0.9338*x4 + x7
     - x9 <= 0.3571;

e3: 0.1238*x1 - 0.004731*x1*x3 + 0.3578*x2*x3 + 0.001637*x2 + 0.9338*x4 - x7
     - x9 <= -0.3571;

e4: 0.2238*x1*x3 + 0.2638*x1 + 0.7623*x2*x3 - 0.07745*x2 - 0.6734*x4 - x7 - x9
     <= 0.6022;

e5: (-0.2238*x1*x3) - 0.2638*x1 - 0.7623*x2*x3 + 0.07745*x2 + 0.6734*x4 + x7
     - x9 <= -0.6022;

e6: x6*x8 + 0.3578*x1 + 0.004731*x2 - x9 <= 0;

e7:  - x6*x8 - 0.3578*x1 - 0.004731*x2 - x9 <= 0;

e8:  - 0.7623*x1 + 0.2238*x2 = -0.3461;

e9: x1^2 + x2^2 - x9 <= 1;

e10: (-x1^2) - x2^2 - x9 <= -1;

e11: x3^2 + x4^2 - x9 <= 1;

e12: (-x3^2) - x4^2 - x9 <= -1;

e13: x5^2 + x6^2 - x9 <= 1;

e14: (-x5^2) - x6^2 - x9 <= -1;

e15: x7^2 + x8^2 - x9 <= 1;

e16: (-x7^2) - x8^2 - x9 <= -1;
