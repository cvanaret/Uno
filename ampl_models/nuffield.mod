#   A problem from economics.
#   Maximize a 2-D integral representing consumer surplus subject to 
#   linear and quadratic constraints representing incentive compatibility

#   Let v( . , . ) : R^2 -> R, Omega = [a,a+1] x [a,a+1], and
#   the corners A, B, C, D be as follows:

#           (a+1,a+1)
#       A *-----* B
#         |     |
#         |     |
#       D *-----* C
#       (a,a)  

#   The problem is to maximize

#      (a+1) line integral_{AB U BC} v(w)dw 
#       - a line integral_{CD U DA} v(w)dw
#       - 3 volume integral_{Omega} v(w)dw

#   subject to v being symmetric (i.e., v(x,y) = v(y,x))
#              v(a,a) = 0
#              nabla_w v(w) >= 0
#              < e, nabla_w v(w) > <= 1
#        and   nabla_ww v(w) positive definite

#   this last constraint is guaranteed by ensuring that

#              d^2 v/dx^2 >= 0
#              d^2 v/dy^2 >= 0
#              ( d^2 v/dx^2 )( d^2 v/dy^2 ) >= ( d^2 v/dxdy )^2

#   Symmetry is ensured by only considering v(x,y) for x <= y

#   Here v(x,y) is the consumer surplus. that is if the consumer values good 
#   1 at x pounds and good 2 at y pounds then they will have a utility 
#   equivalent to v(x,y) pounds after being faced with the optimal monopoly 
#   pricing strategy. (Apparently, from this we can infer what the optimal 
#   pricing strategy was... ).

#   More background is available from

#   "Optimal Selling Strategies: When to haggle, when to hold firm",
#     Riley and Zeckhauser. The Quarterly Journal of Economics, 1983, and

#   "Multidimensional Incentive Compatibility and Mechanism Design", 
#     McAfee and McMillan. The Journal of Economic Theory, 1988.

#   Source: John Thanassoulis <john.thanassoulis@nuffield.oxford.ac.uk>

#   AMPL model: Nick Gould

#   Standard finite-differences are used to ap[proximate derivatives, and 
#   1- and 2-D trapezoidal rules to approximate integrals

param a := 5;  #the bottom left hand point of the grid
param N := 40; #the number of grid points along the bottom row (starting from 0)
param h := 1/N; #the distance between grid points

var v_sym{i in 0..N, j in 0..i}; 

var v {i in 0..N, j in 0..N} = if (j<=i) then v_sym[i,j] else v_sym[j,i];

##dv/dx 

var dvdx {i in 0..N-1, j in 0..N} = (v[i+1,j]-v[i,j])/h;

##dv/dy 

var dvdy {i in 0..N, j in 0..N-1} = (v[i,j+1]-v[i,j])/h;

##d2v/dx2

var d2vdx2 {i in 1..N-1, j in 0..N} = (v[i+1,j]-2*v[i,j]+v[i-1,j])/h^2;

##d2v/dy2

var d2vdy2 {i in 0..N, j in 1..N-1} = (v[i,j+1]-2*v[i,j]+v[i,j-1])/h^2;

##d2v/dxdy

var d2vdxdy {i in 0..N-1, j in 0..N-1} =
    (v[i+1,j+1]-v[i,j+1]-v[i+1,j]+v[i,j])/h^2;

maximize integral: 
    (a+1)*
    h*(sum{j in 1..N-1} (v[N,j]+v[j,N]) + (1/2)*(v[N,0]+2*v[N,N]+v[0,N])) 
    -
    a*
    h*(sum{j in 1..N-1} (v[0,j]+v[j,0]) + (h/2)*(v[N,0]+2*v[0,0]+v[0,N])) 
    -
    3*
    h^2*(
       sum{i in 1..N-1, j in 1..N-1} v[i,j] +
       sum{j in 1..N-1} (v[j,0] + v[j,N] + v[0,j] + v[N,j])/2 +
       (v[0,0] + v[N,0] + v[0,N] + v[N,N])/4
    );

#subject to symmetry{i in 0..N, j in 0..i-1}: v[i,j] = v[j,i];

subject to dx{i in 0..N-1, j in 0..N}:
	dvdx[i,j] >= 0;
	
subject to dy{i in 0..N, j in 0..N-1}:
	dvdy[i,j] >= 0;

subject to prob{i in 0..N-1, j in 0..N-1}:
	(dvdx[i,j]+dvdx[i,j+1])/2 + (dvdy[i,j]+dvdy[i+1,j])/2 <= 1;

subject to d2x{i in 1..N-1, j in 0..N}:
	d2vdx2[i,j] >= 
	1e-8; 
	#0;

#subject to d2y{i in 0..N, j in 1..N-1}:
#	d2vdy2[i,j] >=0;

#subject to conv{i in 1..N-1, j in 1..N-1}:
#	d2vdx2[i,j]*d2vdy2[i,j]-d2vdxdy[i,j]^2 >= 
#	#0; 
#	1.0e-8;

subject to conv{i in 1..N-1, j in 1..N-1}:
	(d2vdy2[i,j])-(d2vdxdy[i,j])^2/d2vdx2[i,j] >= 1e-8; #0; #1.0e-8;

subject to bound: v[0,0]=0;

let {i in 0..N, j in 0..i} v_sym[i,j] := (i^2+j^2)/(4*N^2);

#let {i in 0..N, j in 0..i: i+j>N*0.0985434 && j >  i - N/3} 
#		v_sym[i,j] := (i/N + j/N - 0.0985434)/2;
#let {i in 0..N, j in 0..i: i+j>N*0.0985434 && j <= i - N/3} 
#		v_sym[i,j] := (i/N - 1/3 + (1/3 - 0.0985434)/2);
#let {i in 0..N, j in 0..i: i+j<=N*0.0985434}
#		v_sym[i,j] := 0;
#display integral/2;

solve;

printf {i in 0..N, j in 0..N}: 
    "%10.4f %10.4f %10.4f \n", a+i/N, a+j/N, v[i,j];
