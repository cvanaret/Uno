#  NLP written by GAMS Convert at 06/20/02 11:43:32
#  
#  Equation counts
#     Total       E       G       L       N       X
#         7       1       3       3       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         7       7       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        19      11       8       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 5, >= 0;
var x2 := 1, >= 0;
var x3 := 5, >= 1, <= 5;
var x4 >= 0, <= 6;
var x5 := 5, >= 1, <= 5;
var x6 := 10, >= 0, <= 10;

minimize obj: (-25*(x1 - 2)^2) - (x2 - 2)^2 - (x3 - 1)^2 - (x4 - 4)^2 - (x5 - 1
              )^2 - (x6 - 4)^2;

subject to

e2: (x3 - 3)^2 + x4 >= 4;

e3: (x5 - 3)^2 + x6 >= 4;

e4:    x1 - 3*x2 <= 2;

e5:  - x1 + x2 <= 2;

e6:    x1 + x2 <= 6;

e7:    x1 + x2 >= 2;
