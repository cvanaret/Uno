#  NLP written by GAMS Convert at 07/08/02 14:14:39
#  
#  Equation counts
#     Total       E       G       L       N       X       C
#         5       5       0       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         7       7       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        21       9      12       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x2 := 0.333333333333333;
var x3 := 0.333333333333333;
var x4 := 0.333333333333333;
var x5 := 2;
var x6 := 1;
var x7 := 1;

minimize obj: x2*log(x2) + x3*log(x3) + x4*log(x4) - log(x5 - x7) + x5 - 
              0.353553390593274*x6*log((x5 + 2.41421356237309*x7)/(x5 - 
              0.414213562373095*x7))/x7 + 1.42876598488588*x2 + 
              1.27098480432594*x3 + 1.62663700075562*x4 - 1;

subject to

e2: x5^3 - (1 - x7)*x5^2 + (-3*x7^2 - 2*x7 + x6)*x5 - x6*x7 + x7^3 + x7^2 = 0;

e3:  - (0.142724*x2*x2 + 0.206577*x2*x3 + 0.342119*x2*x4 + 0.206577*x3*x2 + 
    0.323084*x3*x3 + 0.547748*x3*x4 + 0.342119*x4*x2 + 0.547748*x4*x3 + 
    0.968906*x4*x4) + x6 = 0;

e4:  - 0.0815247*x2 - 0.0907391*x3 - 0.13705*x4 + x7 = 0;

e5:    x2 + x3 + x4 = 1;
