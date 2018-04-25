#  NLP written by GAMS Convert at 07/09/02 21:32:09
#  
#  Equation counts
#     Total       E       G       L       N       X       C
#         8       1       0       7       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         5       5       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        17      14       3       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1;
var x2;
var x3;
var x4;

minimize obj:    x4;

subject to

e2: x1^4*x2^4 - x1^4 - x2^4*x3 <= 0;

e3:  - x1 - 0.25*x4 <= -1.4;

e4:    x1 - 0.25*x4 <= 1.4;

e5:  - x2 - 0.2*x4 <= -1.5;

e6:    x2 - 0.2*x4 <= 1.5;

e7:  - x3 - 0.2*x4 <= -0.8;

e8:    x3 - 0.2*x4 <= 0.8;
