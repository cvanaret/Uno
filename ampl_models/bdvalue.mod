# AMPL Model by Hande Y. Benson
#
# Copyright (C) 2001 Princeton University
# All Rights Reserved
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby
# granted, provided that the above copyright notice appear in all
# copies and that the copyright notice and this
# permission notice appear in all supporting documentation.                     

#   Source:  problem 28 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   SIF input: Ph. Toint, Dec 1989.

#   classification NOR2-MN-V-V

param ndp:=5002;
param h:=1/(ndp-1);
var x{i in 1..ndp} := ( (i-1)*h )*( (i-1)*h-1 );

minimize f: 0;
subject to cons{i in 2..ndp-1}: ( -x[i-1]+2*x[i]-x[i+1]+0.5*h^2*(x[i]+i*h+1)^3 ) = 0;

fix x[1] := 0.0;
fix x[ndp] := 0.0;

solve; display f; display x;
