#  NLP written by GAMS Convert at 06/20/02 11:41:12
#  
#  Equation counts
#     Total       E       G       L       N       X
#        12       1       0      11       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        11      11       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#       112     105       7       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 1, >= 0, <= 1;
var x2 := 0.90755, >= 0, <= 1;
var x3 >= 0, <= 1;
var x4 := 1, >= 0, <= 1;
var x5 := 0.71509, >= 0, <= 1;
var x6 := 1, >= 0, <= 1;
var x7 >= 0, <= 1;
var x8 := 0.91698, >= 0, <= 1;
var x9 := 1, >= 0, <= 1;
var x10 := 1, >= 0, <= 1;

minimize obj: (-0.5*(10*x1*x1 + 10*x2*x2 + 10*x3*x3 + 10*x4*x4 + 10*x5*x5 + 10*
              x6*x6 + 10*x7*x7)) - 20*x1 - 80*x2 - 20*x3 - 50*x4 - 60*x5 - 90*
              x6 + 10*x8 + 10*x9 + 10*x10;

subject to

e2:  - 2*x1 - 6*x2 - x3 - 3*x5 - 3*x6 - 2*x7 - 6*x8 - 2*x9 - 2*x10 <= -4;

e3:    6*x1 - 5*x2 + 8*x3 - 3*x4 + x6 + 3*x7 + 8*x8 + 9*x9 - 3*x10 <= 22;

e4:  - 5*x1 + 6*x2 + 5*x3 + 3*x4 + 8*x5 - 8*x6 + 9*x7 + 2*x8 - 9*x10 <= -6;

e5:    9*x1 + 5*x2 - 9*x4 + x5 - 8*x6 + 3*x7 - 9*x8 - 9*x9 - 3*x10 <= -23;

e6:  - 8*x1 + 7*x2 - 4*x3 - 5*x4 - 9*x5 + x6 - 7*x7 - x8 + 3*x9 - 2*x10 <= -12;

e7:  - 7*x1 - 5*x2 - 2*x3 - 6*x5 - 6*x6 - 7*x7 - 6*x8 + 7*x9 + 7*x10 <= -3;

e8:    x1 - 3*x2 - 3*x3 - 4*x4 - x5 - 4*x7 + x8 + 6*x9 <= 1;

e9:    x1 - 2*x2 + 6*x3 + 9*x4 - 7*x6 + 9*x7 - 9*x8 - 6*x9 + 4*x10 <= 12;

e10:  - 4*x1 + 6*x2 + 7*x3 + 2*x4 + 2*x5 + 6*x7 + 6*x8 - 7*x9 + 4*x10 <= 15;

e11:    x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 <= 9;

e12:  - x1 - x2 - x3 - x4 - x5 - x6 - x7 - x8 - x9 - x10 <= -1;
