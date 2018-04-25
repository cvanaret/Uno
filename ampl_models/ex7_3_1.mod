#  NLP written by GAMS Convert at 07/09/02 21:32:08
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


var x1 >= 0;
var x2 >= 0;
var x3 >= 0;
var x4 >= 0;

minimize obj:    x4;

subject to

e2: 10*x2^2*x3^3 + 10*x2^3*x3^2 + 200*x2^2*x3^2 + 100*x2^3*x3 + 100*x2*x3^3 + 
    x1*x2*x3^2 + x1*x2^2*x3 + 1000*x2*x3^2 + 8*x1*x3^2 + 1000*x2^2*x3 + 8*x1*x2
    ^2 + 6*x1*x2*x3 - x1^2 + 60*x1*x3 + 60*x1*x2 - 200*x1 <= 0;

e3:  - x1 - 800*x4 <= -800;

e4:    x1 - 800*x4 <= 800;

e5:  - x2 - 2*x4 <= -4;

e6:    x2 - 2*x4 <= 4;

e7:  - x3 - 3*x4 <= -6;

e8:    x3 - 3*x4 <= 6;
