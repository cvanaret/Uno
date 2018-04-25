#  NLP written by GAMS Convert at 06/20/02 11:59:02
#  
#  Equation counts
#     Total       E       G       L       N       X
#         6       2       0       4       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         6       6       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        20       7      13       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 >= -5, <= 5;
var x2 >= -5, <= 5;
var x3 >= -5, <= 5;
var x4 >= -5, <= 5;
var x5 >= -5, <= 5;

minimize obj: (x1 - 1)^2 + (x1 - x2)^2 + (x2 - x3)^3 + (x3 - x4)^4 + (x4 - x5)^
              4;

subject to

e1: x2^2 + x3^3 + x1 <= 6.24264068711929;

e2: (-x3^3) - x2^2 - x1 <= -6.24264068711929;

e3:  - x3^2 + x2 + x4 <= 0.82842712474619;

e4: x3^2 - x2 - x4 <= -0.82842712474619;

e5: 0.5*x1*x5 + 0.5*x1*x5 = 2;
