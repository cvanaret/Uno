#  NLP written by GAMS Convert at 06/20/02 11:20:18
#  
#  Equation counts
#     Total       E       G       L       N       X
#         8       8       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        15      15       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        32      13      19       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x2 := 1.745, >= 0, <= 2;
var x3 := 1.2, >= 0, <= 1.6;
var x4 := 1.1, >= 0, <= 1.2;
var x5 := 3.048, >= 0, <= 5;
var x6 := 1.974, >= 0, <= 2;
var x7 := 0.893, >= 0.85, <= 0.93;
var x8 := 0.928, >= 0.9, <= 0.95;
var x9 := 8, >= 3, <= 12;
var x10 := 3.6, >= 1.2, <= 4;
var x11 := 1.45, >= 1.45, <= 1.62;
var x12 := 1, >= 0.99, <= 1.01010101010101;
var x13 := 1, >= 0.99, <= 1.01010101010101;
var x14 := 1, >= 0.9, <= 1.11111111111111;
var x15 := 1, >= 0.99, <= 1.01010101010101;

minimize obj:  - 6.3*x5*x8 + 5.04*x2 + 0.35*x3 + x4 + 3.36*x6;

subject to

e2:  - 0.819672131147541*x2 + x5 - 0.819672131147541*x6 = 0;

e3: 0.98*x4 - x7*(0.01*x5*x10 + x4) = 0;

e4:  - x2*x9 + 10*x3 + x6 = 0;

e5: x5*x12 - x2*(1.12 + 0.13167*x9 - 0.0067*x9*x9) = 0;

e6: x8*x13 - 0.01*(1.098*x9 - 0.038*x9*x9) - 0.325*x7 = 0.57425;

e7: x10*x14 + 22.2*x11 = 35.82;

e8: x11*x15 - 3*x8 = -1.33;
