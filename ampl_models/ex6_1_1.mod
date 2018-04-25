#  NLP written by GAMS Convert at 06/20/02 11:51:07
#  
#  Equation counts
#     Total       E       G       L       N       X
#         7       7       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         9       9       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        25       5      20       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x2 := 0.4993, >= 1E-7, <= 0.5;
var x3 := 0.0007, >= 1E-7, <= 0.5;
var x4 := 0.3441, >= 1E-7, <= 0.5;
var x5 := 0.1559, >= 1E-7, <= 0.5;
var x6 := 0.901221308981222, >= 0;
var x7 := 0.0274569351394739, >= 0;
var x8 := 0.691165161172019, >= 0;
var x9 := 0.998619236157215, >= 0;

minimize obj: x2*(log(x2) - log(x2 + x4)) + x4*(log(x4) - log(x2 + x4)) + x3*(
              log(x3) - log(x3 + x5)) + x5*(log(x5) - log(x3 + x5)) + 
              0.925356626778358*x2*x8 + 0.746014540096753*x4*x6 + 
              0.925356626778358*x3*x9 + 0.746014540096753*x5*x7;

subject to

e2: x6*(x2 + 0.159040857374844*x4) - x2 = 0;

e3: x7*(x3 + 0.159040857374844*x5) - x3 = 0;

e4: x8*(0.307941026821595*x2 + x4) - x4 = 0;

e5: x9*(0.307941026821595*x3 + x5) - x5 = 0;

e6:    x2 + x3 = 0.5;

e7:    x4 + x5 = 0.5;
