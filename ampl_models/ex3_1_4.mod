#  NLP written by GAMS Convert at 06/20/02 11:44:02
#  
#  Equation counts
#     Total       E       G       L       N       X
#         4       1       1       2       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         4       4       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        12       9       3       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 0.5, >= 0, <= 2;
var x2 >= 0;
var x3 := 3, >= 0, <= 3;

minimize obj:  - 2*x1 + x2 - x3;

subject to

e2: x1*(4*x1 - 2*x2 + 2*x3) + x2*(2*x2 - 2*x1 - x3) + x3*(2*x1 - x2 + 2*x3) - 
    20*x1 + 9*x2 - 13*x3 >= -24;

e3:    x1 + x2 + x3 <= 4;

e4:    3*x2 + x3 <= 6;
