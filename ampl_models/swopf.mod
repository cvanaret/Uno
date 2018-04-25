#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem:
#   ********
#   An optimal electrical powerflow system design problem from Switzerland.
#   Source:
#   a contribution to fullfill the LANCELOT academic licence agreement.
#   SIF input: R. Bacher, Dept of Electrical Engineering, ETH Zurich, 
#              November 1994.
#   classification LQR2-RN-83-92
#   Number of nodes       =   7
#   Number of branches    =   7
#   Number of generators  =   3
#   Bus oriented data
#   Branch oriented data (including transformers)
#   Generator oriented data
#   Individual generator oriented data
# Bus oriented data
#   Branch oriented data (including transformers)
# Generator oriented data
# Individual generator oriented data
#   Bus oriented data
#   Branch oriented data (including transformers)
#   Bus oriented data
#   Generator oriented data
#   Bus oriented data
#   Branch oriented data (including transformers)
#   Generator oriented data
#   Individual generator oriented data
#    Bus oriented data
#   Branch oriented data (including transformers)
#   Bus oriented data
#   Branch oriented data (including transformers)
#   Solutions
#LO SOLTN               6.78605619D-2
param first := 1;
param nobranches := 7;
param noshunts := 0;
param notrafos := 3;
param nobusses := 7;
param nogen := 3;
param nogenbk := 3;
param noareas := 0;

var ve0001 := 1.0;
var vf0001;
var v20001 >= 0.81 ,  <= 1.21 ,  := 1.0;
var ve0002 := 1.001;
var vf0002;
var v20002 >= 0.81 ,  <= 1.21 ,  := 1.002;
var ve0003 := 1.05;
var vf0003;
var v20003 >= 0.941 ,  <= 1.21 ,  := 1.102;
var ve0004 := 1.001;
var vf0004;
var v20004 >= 0.941 ,  <= 1.21 ,  := 1.002;
var ve0005 := 1.05;
var vf0005;
var v20005 >= 0.941 ,  <= 1.21 ,  := 1.102;
var ve0006 := 1.001;
var vf0006;
var v20006 >= 0.941 ,  <= 1.21 ,  := 1.002;
var ve0007 := 1.001;
var vf0007;
var v20007 >= 0.941 ,  <= 1.21 ,  := 1.002;
var ei0001 := -0.0050;
var fi0001 := 0.066;
var ej0001 := 0.0050;
var fj0001 := -0.066;
var pi0001 := -0.0050;
var qi0001 := -0.066;
var pj0001 := 0.0050;
var qj0001 := 0.066;
var ei0002;
var fi0002 := 0.136;
var ej0002;
var fj0002 := 0.136;
var pi0002;
var qi0002 := -0.136;
var pj0002;
var qj0002 := -0.136;
var ei0003;
var fi0003 := 0.091;
var ej0003;
var fj0003 := 0.091;
var pi0003;
var qi0003 := -0.091;
var pj0003;
var qj0003 := -0.091;
var ei0004 := 0.338;
var fi0004 := -4.055;
var ej0004 := -0.338;
var fj0004 := 4.055;
var pi0004 := 0.355;
var qi0004 := 4.258;
var pj0004 := -0.338;
var qj0004 := -4.059;
var ei0005;
var fi0005 := 0.136;
var ej0005;
var fj0005 := 0.136;
var pi0005;
var qi0005 := -0.136;
var pj0005;
var qj0005 := -0.136;
var ei0006 := 0.169;
var fi0006 := -2.028;
var ej0006 := -0.169;
var fj0006 := 2.028;
var pi0006 := 0.177;
var qi0006 := 2.129;
var pj0006 := -0.169;
var qj0006 := -2.03;
var ei0007;
var fi0007 := 0.091;
var ej0007;
var fj0007 := 0.091;
var pi0007;
var qi0007 := -0.091;
var pj0007;
var qj0007 := -0.091;
var pg0001 >= 0.5 ,  <= 10.0 ,  := 3.0;
var pg0002 >= 0.5 ,  <= 10.0 ,  := 5.0;
var pg0003 >= 0.2 ,  <= 4.0 ,  := 2.0;
var qg0001;
var qg0002;
var qg0003;

minimize obj:
	pi0001 + pj0001 + pi0002 + pj0002 + pi0003 + pj0003 + pi0004 + pj0004 + pi0005 + pj0005 + pi0006 + pj0006 + pi0007 + pj0007;

subject to gv20001:
	-ve0001 * ve0001 - vf0001 * vf0001 + v20001 = 0;
subject to slf0000:
	vf0001 = 0;
subject to gv20002:
	-ve0002 * ve0002 - vf0002 * vf0002 + v20002 = 0;
subject to gv20003:
	-ve0003 * ve0003 - vf0003 * vf0003 + v20003 = 0;
subject to gv20004:
	-ve0004 * ve0004 - vf0004 * vf0004 + v20004 = 0;
subject to gv20005:
	-ve0005 * ve0005 - vf0005 * vf0005 + v20005 = 0;
subject to gv20006:
	-ve0006 * ve0006 - vf0006 * vf0006 + v20006 = 0;
subject to gv20007:
	-ve0007 * ve0007 - vf0007 * vf0007 + v20007 = 0;
subject to gei0001:
	ei0001 - 5.299*ve0001 - 66.243*vf0001 + 5.299*ve0002 + 66.243*vf0002 = 0;
subject to gfi0001:
	fi0001 + 66.243*ve0001 - 5.299*vf0001 - 66.243*ve0002 + 5.299*vf0002 = 0;
subject to gej0001:
	ej0001 - 5.299*ve0002 - 66.243*vf0002 + 5.299*ve0001 + 66.243*vf0001 = 0;
subject to gfj0001:
	fj0001 + 66.243*ve0002 - 5.299*vf0002 - 66.243*ve0001 + 5.299*vf0001 = 0;
subject to gpi0001:
	-ei0001 * ve0001 - fi0001 * vf0001 + pi0001 = 0;
subject to gqi0001:
	-ei0001 * vf0001 + fi0001 * ve0001 + qi0001 = 0;
subject to gpj0001:
	-ej0001 * ve0002 - fj0001 * vf0002 + pj0001 = 0;
subject to gqj0001:
	-ej0001 * vf0002 + fj0001 * ve0002 + qj0001 = 0;
subject to gpni0001:
	-pi0001 + pg0001 = 0;
subject to gpni0002:
	-pj0001 - pi0002 - pi0003 - 2.0 = 0;
subject to gqni0001:
	-qi0001 + qg0001 = 0;
subject to gqni0002:
	-qj0001 - qi0002 - qi0003 - 3.0 = 0;
subject to gmxi0001:
	0 >= pi0001 * pi0001 + qi0001 * qi0001 - 16.0;
subject to gmxj0001:
	0 >= pj0001 * pj0001 + qj0001 * qj0001 - 16.0;
subject to gei0002:
	ei0002 - 1.175*ve0002 - 6.915*vf0002 + 1.175*ve0006 + 7.051*vf0006 = 0;
subject to gfi0002:
	fi0002 + 6.915*ve0002 - 1.175*vf0002 - 7.051*ve0006 + 1.175*vf0006 = 0;
subject to gej0002:
	ej0002 - 1.175*ve0006 - 6.915*vf0006 + 1.175*ve0002 + 7.051*vf0002 = 0;
subject to gfj0002:
	fj0002 + 6.915*ve0006 - 1.175*vf0006 - 7.051*ve0002 + 1.175*vf0002 = 0;
subject to gpi0002:
	-ei0002 * ve0002 - fi0002 * vf0002 + pi0002 = 0;
subject to gqi0002:
	-ei0002 * vf0002 + fi0002 * ve0002 + qi0002 = 0;
subject to gpj0002:
	-ej0002 * ve0006 - fj0002 * vf0006 + pj0002 = 0;
subject to gqj0002:
	-ej0002 * vf0006 + fj0002 * ve0006 + qj0002 = 0;
subject to gpni0006:
	-pj0002 - pj0006 - pi0007 - 1.0 = 0;
subject to gqni0006:
	-qj0002 - qj0006 - qi0007 - 0.3 = 0;
subject to gmxi0002:
	0 >= pi0002 * pi0002 + qi0002 * qi0002 - 4.0;
subject to gmxj0002:
	0 >= pj0002 * pj0002 + qj0002 * qj0002 - 4.0;
subject to gei0003:
	ei0003 - 1.726*ve0002 - 10.498*vf0002 + 1.726*ve0004 + 10.588*vf0004 = 0;
subject to gfi0003:
	fi0003 + 10.498*ve0002 - 1.726*vf0002 - 10.588*ve0004 + 1.726*vf0004 = 0;
subject to gej0003:
	ej0003 - 1.726*ve0004 - 10.498*vf0004 + 1.726*ve0002 + 10.588*vf0002 = 0;
subject to gfj0003:
	fj0003 + 10.498*ve0004 - 1.726*vf0004 - 10.588*ve0002 + 1.726*vf0002 = 0;
subject to gpi0003:
	-ei0003 * ve0002 - fi0003 * vf0002 + pi0003 = 0;
subject to gqi0003:
	-ei0003 * vf0002 + fi0003 * ve0002 + qi0003 = 0;
subject to gpj0003:
	-ej0003 * ve0004 - fj0003 * vf0004 + pj0003 = 0;
subject to gqj0003:
	-ej0003 * vf0004 + fj0003 * ve0004 + qj0003 = 0;
subject to gpni0004:
	-pj0003 - pj0004 - pi0005 - 2.0 = 0;
subject to gqni0004:
	-qj0003 - qj0004 - qi0005 - 0.2 = 0;
subject to gmxi0003:
	0 >= pi0003 * pi0003 + qi0003 * qi0003 - 4.0;
subject to gmxj0003:
	0 >= pj0003 * pj0003 + qj0003 * qj0003 - 4.0;
subject to gei0004:
	ei0004 - 6.897*ve0003 - 82.759*vf0003 + 6.897*ve0004 + 82.759*vf0004 = 0;
subject to gfi0004:
	fi0004 + 82.759*ve0003 - 6.897*vf0003 - 82.759*ve0004 + 6.897*vf0004 = 0;
subject to gej0004:
	ej0004 - 6.897*ve0004 - 82.759*vf0004 + 6.897*ve0003 + 82.759*vf0003 = 0;
subject to gfj0004:
	fj0004 + 82.759*ve0004 - 6.897*vf0004 - 82.759*ve0003 + 6.897*vf0003 = 0;
subject to gpi0004:
	-ei0004 * ve0003 - fi0004 * vf0003 + pi0004 = 0;
subject to gqi0004:
	-ei0004 * vf0003 + fi0004 * ve0003 + qi0004 = 0;
subject to gpj0004:
	-ej0004 * ve0004 - fj0004 * vf0004 + pj0004 = 0;
subject to gqj0004:
	-ej0004 * vf0004 + fj0004 * ve0004 + qj0004 = 0;
subject to gpni0003:
	-pi0004 + pg0002 - 0.6 = 0;
subject to gqni0003:
	-qi0004 + qg0002 - 0.08 = 0;
subject to gmxi0004:
	0 >= pi0004 * pi0004 + qi0004 * qi0004 - 25.0;
subject to gmxj0004:
	0 >= pj0004 * pj0004 + qj0004 * qj0004 - 25.0;
subject to gei0005:
	ei0005 - 1.175*ve0004 - 6.915*vf0004 + 1.175*ve0007 + 7.051*vf0007 = 0;
subject to gfi0005:
	fi0005 + 6.915*ve0004 - 1.175*vf0004 - 7.051*ve0007 + 1.175*vf0007 = 0;
subject to gej0005:
	ej0005 - 1.175*ve0007 - 6.915*vf0007 + 1.175*ve0004 + 7.051*vf0004 = 0;
subject to gfj0005:
	fj0005 + 6.915*ve0007 - 1.175*vf0007 - 7.051*ve0004 + 1.175*vf0004 = 0;
subject to gpi0005:
	-ei0005 * ve0004 - fi0005 * vf0004 + pi0005 = 0;
subject to gqi0005:
	-ei0005 * vf0004 + fi0005 * ve0004 + qi0005 = 0;
subject to gpj0005:
	-ej0005 * ve0007 - fj0005 * vf0007 + pj0005 = 0;
subject to gqj0005:
	-ej0005 * vf0007 + fj0005 * ve0007 + qj0005 = 0;
subject to gpni0007:
	-pj0005 - pj0007 - 2.0 = 0;
subject to gqni0007:
	-qj0005 - qj0007 - 1.0 = 0;
subject to gmxi0005:
	0 >= pi0005 * pi0005 + qi0005 * qi0005 - 4.0;
subject to gmxj0005:
	0 >= pj0005 * pj0005 + qj0005 * qj0005 - 4.0;
subject to gei0006:
	ei0006 - 3.448*ve0005 - 41.379*vf0005 + 3.448*ve0006 + 41.379*vf0006 = 0;
subject to gfi0006:
	fi0006 + 41.379*ve0005 - 3.448*vf0005 - 41.379*ve0006 + 3.448*vf0006 = 0;
subject to gej0006:
	ej0006 - 3.448*ve0006 - 41.379*vf0006 + 3.448*ve0005 + 41.379*vf0005 = 0;
subject to gfj0006:
	fj0006 + 41.379*ve0006 - 3.448*vf0006 - 41.379*ve0005 + 3.448*vf0005 = 0;
subject to gpi0006:
	-ei0006 * ve0005 - fi0006 * vf0005 + pi0006 = 0;
subject to gqi0006:
	-ei0006 * vf0005 + fi0006 * ve0005 + qi0006 = 0;
subject to gpj0006:
	-ej0006 * ve0006 - fj0006 * vf0006 + pj0006 = 0;
subject to gqj0006:
	-ej0006 * vf0006 + fj0006 * ve0006 + qj0006 = 0;
subject to gpni0005:
	-pi0006 + pg0003 - 0.5 = 0;
subject to gqni0005:
	-qi0006 + qg0003 - 0.05 = 0;
subject to gmxi0006:
	0 >= pi0006 * pi0006 + qi0006 * qi0006 - 6.25;
subject to gmxj0006:
	0 >= pj0006 * pj0006 + qj0006 * qj0006 - 6.25;
subject to gei0007:
	ei0007 - 1.726*ve0006 - 10.498*vf0006 + 1.726*ve0007 + 10.588*vf0007 = 0;
subject to gfi0007:
	fi0007 + 10.498*ve0006 - 1.726*vf0006 - 10.588*ve0007 + 1.726*vf0007 = 0;
subject to gej0007:
	ej0007 - 1.726*ve0007 - 10.498*vf0007 + 1.726*ve0006 + 10.588*vf0006 = 0;
subject to gfj0007:
	fj0007 + 10.498*ve0007 - 1.726*vf0007 - 10.588*ve0006 + 1.726*vf0006 = 0;
subject to gpi0007:
	-ei0007 * ve0006 - fi0007 * vf0006 + pi0007 = 0;
subject to gqi0007:
	-ei0007 * vf0006 + fi0007 * ve0006 + qi0007 = 0;
subject to gpj0007:
	-ej0007 * ve0007 - fj0007 * vf0007 + pj0007 = 0;
subject to gqj0007:
	-ej0007 * vf0007 + fj0007 * ve0007 + qj0007 = 0;
subject to gmxi0007:
	0 >= pi0007 * pi0007 + qi0007 * qi0007 - 4.0;
subject to gmxj0007:
	0 >= pj0007 * pj0007 + qj0007 * qj0007 - 4.0;

solve;
display ve0001;
display vf0001;
display v20001;
display ve0002;
display vf0002;
display v20002;
display ve0003;
display vf0003;
display v20003;
display ve0004;
display vf0004;
display v20004;
display ve0005;
display vf0005;
display v20005;
display ve0006;
display vf0006;
display v20006;
display ve0007;
display vf0007;
display v20007;
display ei0001;
display fi0001;
display ej0001;
display fj0001;
display pi0001;
display qi0001;
display pj0001;
display qj0001;
display ei0002;
display fi0002;
display ej0002;
display fj0002;
display pi0002;
display qi0002;
display pj0002;
display qj0002;
display ei0003;
display fi0003;
display ej0003;
display fj0003;
display pi0003;
display qi0003;
display pj0003;
display qj0003;
display ei0004;
display fi0004;
display ej0004;
display fj0004;
display pi0004;
display qi0004;
display pj0004;
display qj0004;
display ei0005;
display fi0005;
display ej0005;
display fj0005;
display pi0005;
display qi0005;
display pj0005;
display qj0005;
display ei0006;
display fi0006;
display ej0006;
display fj0006;
display pi0006;
display qi0006;
display pj0006;
display qj0006;
display ei0007;
display fi0007;
display ej0007;
display fj0007;
display pi0007;
display qi0007;
display pj0007;
display qj0007;
display pg0001;
display pg0002;
display pg0003;
display qg0001;
display qg0002;
display qg0003;
display obj;
