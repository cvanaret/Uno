param N := 10;
var x{i in 1..N} >= (-1)^i*0.2*i, <= 0.8*i, := 1.;

minimize f:
	sum{i in 1..N-1} (100*(x[i + 1] - x[i]^2)^2 + (1 - x[i])^2);

