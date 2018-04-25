#  NLP written by GAMS Convert at 06/20/02 13:05:45
#  
#  Equation counts
#     Total       E       G       L       N       X
#         6       6       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         6       6       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        20      10      10       0


var x1 := 1;
var x2 := 1;
var x3 := 1;
var x4 := 1;
var x5 := 1;
var x6 := 1;

minimize obj: x1;

subject to

e1: x1*x2 = 1;

e2: x3/x1/x4 = 4.8;

e3: x5/x2/x6 = 0.98;

e4: x6*x4 = 1;

e5:    x1 - x2 + 1E-7*x3 - 1E-5*x5 = 0;

e6:    2*x1 - 2*x2 + 1E-7*x3 - 0.01*x4 - 1E-5*x5 + 0.01*x6 = 0;
