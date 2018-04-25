#  NLP written by GAMS Convert at 07/09/02 21:32:10
#  
#  Equation counts
#     Total       E       G       L       N       X       C
#         9       3       0       6       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         6       6       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        22      17       5       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1;
var x2;
var x3;
var x4 >= 0, <= 10;
var x5;

minimize obj:    x5;

subject to

e2: 9.625*x1*x4 - 4*x1 - 78*x4 + 16*x2*x4 - x2 + 16*x4^2 + x3 = -12;

e3: 16*x1*x4 - 19*x1 - 24*x4 - 8*x2 - x3 = -44;

e4:    x1 - 0.25*x5 <= 2.25;

e5:  - x1 - 0.25*x5 <= -2.25;

e6:  - x2 - 0.5*x5 <= -1.5;

e7:    x2 - 0.5*x5 <= 1.5;

e8:  - x3 - 1.5*x5 <= -1.5;

e9:    x3 - 1.5*x5 <= 1.5;
