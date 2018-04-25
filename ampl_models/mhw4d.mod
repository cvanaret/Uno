#  NLP written by GAMS Convert at 06/20/02 12:50:59
#  
#  Equation counts
#     Total       E       G       L       N       X
#         4       4       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         6       6       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        14       4      10       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x2 := -1;
var x3 := 2;
var x4 := 1;
var x5 := -2;
var x6 := -2;

minimize obj: (x2 - 1)^2 + (x2 - x3)^2 + (x3 - x4)^3 + (x4 - x5)^4 + (x5 - x6)^
              4;

subject to

e2: x3^2 + x4^3 + x2 = 6.24264068711929;

e3:  - x4^2 + x3 + x5 = 0.82842712474619;

e4: x2*x6 = 2;
