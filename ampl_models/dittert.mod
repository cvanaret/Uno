#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   Minimize the Dittert function.
#   Source: See Minc, Linear and Multilinear Algebra 21, 1987
#   SIF input: N. Gould, March 1992.
#              minor correction by Ph. Shott, Jan 1995.
#   classification OQR2-AN-V-V
#   Size of matrix
#IE N                   2
#IE N                   3
#IE N                   4
#IE N                   5
#IE N                   6
#IE N                   7
#IE N                   9
#IE N                   10
#   Define constants
#   Compute the number of sub-permanents
#   Sub permanents
#   Entries in the coefficient matrix
#   Row and column sums
#   Define objective function groups
#   Linear terms in the sub-permanent constraints.
#   Row and column sum constraints.
#  Constraint that the sum of all the entries is n.
#  Constraint that the sum of all the entries is n.
#   Entries in the coefficient matrix.
#   Set up the elements associated with sub-permanent constraint K
#   Construct the I-th component of the binary representation of K.
#   Associate elements with nonzero entries in the binary string
#   This corresponds to finding the sub-permanents which occur
#   in the usual expansion of the sub-permanent in terms of its
#   sub-sub-permanents.
#   Solution
#LO SOLTN(2)           5.0D-1
#LO SOLTN(3)           2.22222222D-1
#LO SOLTN(4)           9.375-2
#LO SOLTN(5)           3.84D-2
#LO SOLTN(6)           1.54321098D-2
#LO SOLTN(7)           6.11989902D-3
#LO SOLTN(8)           2.40325927D-3
#LO SOLTN(9)           9.36656708D-4
#LO SOLTN(10)          3.6288D-4
#LO SOLTN(11)          1.39905948D-4
#LO SOLTN(12)          5.37232170D-5
#LO SOLTN(13)          2.05596982D-5
#LO SOLTN(14)          7.84541375D-6
#LO SOLTN(15)          2.98628137D-6
#LO SOLTN(16)          1.13422671D-6
#LO SOLTN(17)          4.29968709D-7
#LO SOLTN(18)          1.62718123D-7
#LO SOLTN(19)          6.14859946D-8
param n := 8;
param rn := 8.0;
param np1 := 1 + (8);
param nmip1 := 8;
param im1 := -1 + (8);
param r2ttn := 128.0;
param s8 := 0.1 + (1.0);
param t0 := 0.1 + (1.0);
param s7 := 0.1 + (2.0);
param t1 := 0.1 + (2.0);
param s6 := 0.1 + (4.0);
param t2 := 0.1 + (4.0);
param s5 := 0.1 + (8.0);
param t3 := 0.1 + (8.0);
param s4 := 0.1 + (16.0);
param t4 := 0.1 + (16.0);
param s3 := 0.1 + (32.0);
param t5 := 0.1 + (32.0);
param s2 := 0.1 + (64.0);
param t6 := 0.1 + (64.0);
param s1 := 0.1 + (128.0);
param t7 := 0.1 + (128.0);
param nm1 := -1 + (8);
param rk1 := 0.1 + (128.0);
param k1 := 1 + (round(0.1 + (128.0)));
param k2 := -1 + (2 * (round(0.1 + (128.0))));
param id := round(0.1 + (8.0));
param pt := round(0.1 + (1.0));
param kk := round(0.1 + (1.0));
param si := 0.1 + (1.0);
param isi := round(0.1 + (1.0));
param bi := round(0.1 + (1.0));
param bisi := round(0.1 + (1.0));
param ri := 8;
param i1 := round(0.1 + (1.0));
param i2 := 0;
param idm2 := -2 + round(0.1 + (1.0));
param rj := 0.1 + (8);
param ipp := round(0.1 + (1.0));
param rjj := 0.1 + (2);
param jj := round(0.1 + (2));
param rd := 0.1 + (8.0);
param rnz1 := 0.1 + (1);
param rnz2 := 0.1 + (2);
param d3 := 0.1 + (2.0);
param d5 := 0.1 + (2.0);
param rnz3 := 0.1 + (3);
param d6 := 0.1 + (2.0);
param d7 := 0.1 + (3.0);
param d9 := 0.1 + (2.0);
param d10 := 0.1 + (2.0);
param d11 := 0.1 + (3.0);
param d12 := 0.1 + (2.0);
param d13 := 0.1 + (3.0);
param rnz4 := 0.1 + (4);
param d14 := 0.1 + (3.0);
param d15 := 0.1 + (4.0);
param d17 := 0.1 + (2.0);
param d18 := 0.1 + (2.0);
param d19 := 0.1 + (3.0);
param d20 := 0.1 + (2.0);
param d21 := 0.1 + (3.0);
param d22 := 0.1 + (3.0);
param d23 := 0.1 + (4.0);
param d24 := 0.1 + (2.0);
param d25 := 0.1 + (3.0);
param d26 := 0.1 + (3.0);
param d27 := 0.1 + (4.0);
param d28 := 0.1 + (3.0);
param d29 := 0.1 + (4.0);
param rnz5 := 0.1 + (5);
param d30 := 0.1 + (4.0);
param d31 := 0.1 + (5.0);
param d33 := 0.1 + (2.0);
param d34 := 0.1 + (2.0);
param d35 := 0.1 + (3.0);
param d36 := 0.1 + (2.0);
param d37 := 0.1 + (3.0);
param d38 := 0.1 + (3.0);
param d39 := 0.1 + (4.0);
param d40 := 0.1 + (2.0);
param d41 := 0.1 + (3.0);
param d42 := 0.1 + (3.0);
param d43 := 0.1 + (4.0);
param d44 := 0.1 + (3.0);
param d45 := 0.1 + (4.0);
param d46 := 0.1 + (4.0);
param d47 := 0.1 + (5.0);
param d48 := 0.1 + (2.0);
param d49 := 0.1 + (3.0);
param d50 := 0.1 + (3.0);
param d51 := 0.1 + (4.0);
param d52 := 0.1 + (3.0);
param d53 := 0.1 + (4.0);
param d54 := 0.1 + (4.0);
param d55 := 0.1 + (5.0);
param d56 := 0.1 + (3.0);
param d57 := 0.1 + (4.0);
param d58 := 0.1 + (4.0);
param d59 := 0.1 + (5.0);
param d60 := 0.1 + (4.0);
param d61 := 0.1 + (5.0);
param rnz6 := 0.1 + (6);
param d62 := 0.1 + (5.0);
param d63 := 0.1 + (6.0);
param d65 := 0.1 + (2.0);
param d66 := 0.1 + (2.0);
param d67 := 0.1 + (3.0);
param d68 := 0.1 + (2.0);
param d69 := 0.1 + (3.0);
param d70 := 0.1 + (3.0);
param d71 := 0.1 + (4.0);
param d72 := 0.1 + (2.0);
param d73 := 0.1 + (3.0);
param d74 := 0.1 + (3.0);
param d75 := 0.1 + (4.0);
param d76 := 0.1 + (3.0);
param d77 := 0.1 + (4.0);
param d78 := 0.1 + (4.0);
param d79 := 0.1 + (5.0);
param d80 := 0.1 + (2.0);
param d81 := 0.1 + (3.0);
param d82 := 0.1 + (3.0);
param d83 := 0.1 + (4.0);
param d84 := 0.1 + (3.0);
param d85 := 0.1 + (4.0);
param d86 := 0.1 + (4.0);
param d87 := 0.1 + (5.0);
param d88 := 0.1 + (3.0);
param d89 := 0.1 + (4.0);
param d90 := 0.1 + (4.0);
param d91 := 0.1 + (5.0);
param d92 := 0.1 + (4.0);
param d93 := 0.1 + (5.0);
param d94 := 0.1 + (5.0);
param d95 := 0.1 + (6.0);
param d96 := 0.1 + (2.0);
param d97 := 0.1 + (3.0);
param d98 := 0.1 + (3.0);
param d99 := 0.1 + (4.0);
param d100 := 0.1 + (3.0);
param d101 := 0.1 + (4.0);
param d102 := 0.1 + (4.0);
param d103 := 0.1 + (5.0);
param d104 := 0.1 + (3.0);
param d105 := 0.1 + (4.0);
param d106 := 0.1 + (4.0);
param d107 := 0.1 + (5.0);
param d108 := 0.1 + (4.0);
param d109 := 0.1 + (5.0);
param d110 := 0.1 + (5.0);
param d111 := 0.1 + (6.0);
param d112 := 0.1 + (3.0);
param d113 := 0.1 + (4.0);
param d114 := 0.1 + (4.0);
param d115 := 0.1 + (5.0);
param d116 := 0.1 + (4.0);
param d117 := 0.1 + (5.0);
param d118 := 0.1 + (5.0);
param d119 := 0.1 + (6.0);
param d120 := 0.1 + (4.0);
param d121 := 0.1 + (5.0);
param d122 := 0.1 + (5.0);
param d123 := 0.1 + (6.0);
param d124 := 0.1 + (5.0);
param d125 := 0.1 + (6.0);
param rnz7 := 0.1 + (7);
param d126 := 0.1 + (6.0);
param d127 := 0.1 + (7.0);
param d129 := 0.1 + (2.0);
param d130 := 0.1 + (2.0);
param d131 := 0.1 + (3.0);
param d132 := 0.1 + (2.0);
param d133 := 0.1 + (3.0);
param d134 := 0.1 + (3.0);
param d135 := 0.1 + (4.0);
param d136 := 0.1 + (2.0);
param d137 := 0.1 + (3.0);
param d138 := 0.1 + (3.0);
param d139 := 0.1 + (4.0);
param d140 := 0.1 + (3.0);
param d141 := 0.1 + (4.0);
param d142 := 0.1 + (4.0);
param d143 := 0.1 + (5.0);
param d144 := 0.1 + (2.0);
param d145 := 0.1 + (3.0);
param d146 := 0.1 + (3.0);
param d147 := 0.1 + (4.0);
param d148 := 0.1 + (3.0);
param d149 := 0.1 + (4.0);
param d150 := 0.1 + (4.0);
param d151 := 0.1 + (5.0);
param d152 := 0.1 + (3.0);
param d153 := 0.1 + (4.0);
param d154 := 0.1 + (4.0);
param d155 := 0.1 + (5.0);
param d156 := 0.1 + (4.0);
param d157 := 0.1 + (5.0);
param d158 := 0.1 + (5.0);
param d159 := 0.1 + (6.0);
param d160 := 0.1 + (2.0);
param d161 := 0.1 + (3.0);
param d162 := 0.1 + (3.0);
param d163 := 0.1 + (4.0);
param d164 := 0.1 + (3.0);
param d165 := 0.1 + (4.0);
param d166 := 0.1 + (4.0);
param d167 := 0.1 + (5.0);
param d168 := 0.1 + (3.0);
param d169 := 0.1 + (4.0);
param d170 := 0.1 + (4.0);
param d171 := 0.1 + (5.0);
param d172 := 0.1 + (4.0);
param d173 := 0.1 + (5.0);
param d174 := 0.1 + (5.0);
param d175 := 0.1 + (6.0);
param d176 := 0.1 + (3.0);
param d177 := 0.1 + (4.0);
param d178 := 0.1 + (4.0);
param d179 := 0.1 + (5.0);
param d180 := 0.1 + (4.0);
param d181 := 0.1 + (5.0);
param d182 := 0.1 + (5.0);
param d183 := 0.1 + (6.0);
param d184 := 0.1 + (4.0);
param d185 := 0.1 + (5.0);
param d186 := 0.1 + (5.0);
param d187 := 0.1 + (6.0);
param d188 := 0.1 + (5.0);
param d189 := 0.1 + (6.0);
param d190 := 0.1 + (6.0);
param d191 := 0.1 + (7.0);
param d192 := 0.1 + (2.0);
param d193 := 0.1 + (3.0);
param d194 := 0.1 + (3.0);
param d195 := 0.1 + (4.0);
param d196 := 0.1 + (3.0);
param d197 := 0.1 + (4.0);
param d198 := 0.1 + (4.0);
param d199 := 0.1 + (5.0);
param d200 := 0.1 + (3.0);
param d201 := 0.1 + (4.0);
param d202 := 0.1 + (4.0);
param d203 := 0.1 + (5.0);
param d204 := 0.1 + (4.0);
param d205 := 0.1 + (5.0);
param d206 := 0.1 + (5.0);
param d207 := 0.1 + (6.0);
param d208 := 0.1 + (3.0);
param d209 := 0.1 + (4.0);
param d210 := 0.1 + (4.0);
param d211 := 0.1 + (5.0);
param d212 := 0.1 + (4.0);
param d213 := 0.1 + (5.0);
param d214 := 0.1 + (5.0);
param d215 := 0.1 + (6.0);
param d216 := 0.1 + (4.0);
param d217 := 0.1 + (5.0);
param d218 := 0.1 + (5.0);
param d219 := 0.1 + (6.0);
param d220 := 0.1 + (5.0);
param d221 := 0.1 + (6.0);
param d222 := 0.1 + (6.0);
param d223 := 0.1 + (7.0);
param d224 := 0.1 + (3.0);
param d225 := 0.1 + (4.0);
param d226 := 0.1 + (4.0);
param d227 := 0.1 + (5.0);
param d228 := 0.1 + (4.0);
param d229 := 0.1 + (5.0);
param d230 := 0.1 + (5.0);
param d231 := 0.1 + (6.0);
param d232 := 0.1 + (4.0);
param d233 := 0.1 + (5.0);
param d234 := 0.1 + (5.0);
param d235 := 0.1 + (6.0);
param d236 := 0.1 + (5.0);
param d237 := 0.1 + (6.0);
param d238 := 0.1 + (6.0);
param d239 := 0.1 + (7.0);
param d240 := 0.1 + (4.0);
param d241 := 0.1 + (5.0);
param d242 := 0.1 + (5.0);
param d243 := 0.1 + (6.0);
param d244 := 0.1 + (5.0);
param d245 := 0.1 + (6.0);
param d246 := 0.1 + (6.0);
param d247 := 0.1 + (7.0);
param d248 := 0.1 + (5.0);
param d249 := 0.1 + (6.0);
param d250 := 0.1 + (6.0);
param d251 := 0.1 + (7.0);
param d252 := 0.1 + (6.0);
param d253 := 0.1 + (7.0);
param rnz8 := 0.1 + (8);
param d254 := 0.1 + (7.0);
param d255 := 0.1 + (8.0);

var p3 >= 0.0, := 0.0;
var p5 >= 0.0, := 0.0;
var p6 >= 0.0, := 0.0;
var p7 >= 0.0, := 0.0;
var p9 >= 0.0, := 0.0;
var p10 >= 0.0, := 0.0;
var p11 >= 0.0, := 0.0;
var p12 >= 0.0, := 0.0;
var p13 >= 0.0, := 0.0;
var p14 >= 0.0, := 0.0;
var p15 >= 0.0, := 0.0;
var p17 >= 0.0, := 0.0;
var p18 >= 0.0, := 0.0;
var p19 >= 0.0, := 0.0;
var p20 >= 0.0, := 0.0;
var p21 >= 0.0, := 0.0;
var p22 >= 0.0, := 0.0;
var p23 >= 0.0, := 0.0;
var p24 >= 0.0, := 0.0;
var p25 >= 0.0, := 0.0;
var p26 >= 0.0, := 0.0;
var p27 >= 0.0, := 0.0;
var p28 >= 0.0, := 0.0;
var p29 >= 0.0, := 0.0;
var p30 >= 0.0, := 0.0;
var p31 >= 0.0, := 0.0;
var p33 >= 0.0, := 0.0;
var p34 >= 0.0, := 0.0;
var p35 >= 0.0, := 0.0;
var p36 >= 0.0, := 0.0;
var p37 >= 0.0, := 0.0;
var p38 >= 0.0, := 0.0;
var p39 >= 0.0, := 0.0;
var p40 >= 0.0, := 0.0;
var p41 >= 0.0, := 0.0;
var p42 >= 0.0, := 0.0;
var p43 >= 0.0, := 0.0;
var p44 >= 0.0, := 0.0;
var p45 >= 0.0, := 0.0;
var p46 >= 0.0, := 0.0;
var p47 >= 0.0, := 0.0;
var p48 >= 0.0, := 0.0;
var p49 >= 0.0, := 0.0;
var p50 >= 0.0, := 0.0;
var p51 >= 0.0, := 0.0;
var p52 >= 0.0, := 0.0;
var p53 >= 0.0, := 0.0;
var p54 >= 0.0, := 0.0;
var p55 >= 0.0, := 0.0;
var p56 >= 0.0, := 0.0;
var p57 >= 0.0, := 0.0;
var p58 >= 0.0, := 0.0;
var p59 >= 0.0, := 0.0;
var p60 >= 0.0, := 0.0;
var p61 >= 0.0, := 0.0;
var p62 >= 0.0, := 0.0;
var p63 >= 0.0, := 0.0;
var p65 >= 0.0, := 0.0;
var p66 >= 0.0, := 0.0;
var p67 >= 0.0, := 0.0;
var p68 >= 0.0, := 0.0;
var p69 >= 0.0, := 0.0;
var p70 >= 0.0, := 0.0;
var p71 >= 0.0, := 0.0;
var p72 >= 0.0, := 0.0;
var p73 >= 0.0, := 0.0;
var p74 >= 0.0, := 0.0;
var p75 >= 0.0, := 0.0;
var p76 >= 0.0, := 0.0;
var p77 >= 0.0, := 0.0;
var p78 >= 0.0, := 0.0;
var p79 >= 0.0, := 0.0;
var p80 >= 0.0, := 0.0;
var p81 >= 0.0, := 0.0;
var p82 >= 0.0, := 0.0;
var p83 >= 0.0, := 0.0;
var p84 >= 0.0, := 0.0;
var p85 >= 0.0, := 0.0;
var p86 >= 0.0, := 0.0;
var p87 >= 0.0, := 0.0;
var p88 >= 0.0, := 0.0;
var p89 >= 0.0, := 0.0;
var p90 >= 0.0, := 0.0;
var p91 >= 0.0, := 0.0;
var p92 >= 0.0, := 0.0;
var p93 >= 0.0, := 0.0;
var p94 >= 0.0, := 0.0;
var p95 >= 0.0, := 0.0;
var p96 >= 0.0, := 0.0;
var p97 >= 0.0, := 0.0;
var p98 >= 0.0, := 0.0;
var p99 >= 0.0, := 0.0;
var p100 >= 0.0, := 0.0;
var p101 >= 0.0, := 0.0;
var p102 >= 0.0, := 0.0;
var p103 >= 0.0, := 0.0;
var p104 >= 0.0, := 0.0;
var p105 >= 0.0, := 0.0;
var p106 >= 0.0, := 0.0;
var p107 >= 0.0, := 0.0;
var p108 >= 0.0, := 0.0;
var p109 >= 0.0, := 0.0;
var p110 >= 0.0, := 0.0;
var p111 >= 0.0, := 0.0;
var p112 >= 0.0, := 0.0;
var p113 >= 0.0, := 0.0;
var p114 >= 0.0, := 0.0;
var p115 >= 0.0, := 0.0;
var p116 >= 0.0, := 0.0;
var p117 >= 0.0, := 0.0;
var p118 >= 0.0, := 0.0;
var p119 >= 0.0, := 0.0;
var p120 >= 0.0, := 0.0;
var p121 >= 0.0, := 0.0;
var p122 >= 0.0, := 0.0;
var p123 >= 0.0, := 0.0;
var p124 >= 0.0, := 0.0;
var p125 >= 0.0, := 0.0;
var p126 >= 0.0, := 0.0;
var p127 >= 0.0, := 0.0;
var p129 >= 0.0, := 0.0;
var p130 >= 0.0, := 0.0;
var p131 >= 0.0, := 0.0;
var p132 >= 0.0, := 0.0;
var p133 >= 0.0, := 0.0;
var p134 >= 0.0, := 0.0;
var p135 >= 0.0, := 0.0;
var p136 >= 0.0, := 0.0;
var p137 >= 0.0, := 0.0;
var p138 >= 0.0, := 0.0;
var p139 >= 0.0, := 0.0;
var p140 >= 0.0, := 0.0;
var p141 >= 0.0, := 0.0;
var p142 >= 0.0, := 0.0;
var p143 >= 0.0, := 0.0;
var p144 >= 0.0, := 0.0;
var p145 >= 0.0, := 0.0;
var p146 >= 0.0, := 0.0;
var p147 >= 0.0, := 0.0;
var p148 >= 0.0, := 0.0;
var p149 >= 0.0, := 0.0;
var p150 >= 0.0, := 0.0;
var p151 >= 0.0, := 0.0;
var p152 >= 0.0, := 0.0;
var p153 >= 0.0, := 0.0;
var p154 >= 0.0, := 0.0;
var p155 >= 0.0, := 0.0;
var p156 >= 0.0, := 0.0;
var p157 >= 0.0, := 0.0;
var p158 >= 0.0, := 0.0;
var p159 >= 0.0, := 0.0;
var p160 >= 0.0, := 0.0;
var p161 >= 0.0, := 0.0;
var p162 >= 0.0, := 0.0;
var p163 >= 0.0, := 0.0;
var p164 >= 0.0, := 0.0;
var p165 >= 0.0, := 0.0;
var p166 >= 0.0, := 0.0;
var p167 >= 0.0, := 0.0;
var p168 >= 0.0, := 0.0;
var p169 >= 0.0, := 0.0;
var p170 >= 0.0, := 0.0;
var p171 >= 0.0, := 0.0;
var p172 >= 0.0, := 0.0;
var p173 >= 0.0, := 0.0;
var p174 >= 0.0, := 0.0;
var p175 >= 0.0, := 0.0;
var p176 >= 0.0, := 0.0;
var p177 >= 0.0, := 0.0;
var p178 >= 0.0, := 0.0;
var p179 >= 0.0, := 0.0;
var p180 >= 0.0, := 0.0;
var p181 >= 0.0, := 0.0;
var p182 >= 0.0, := 0.0;
var p183 >= 0.0, := 0.0;
var p184 >= 0.0, := 0.0;
var p185 >= 0.0, := 0.0;
var p186 >= 0.0, := 0.0;
var p187 >= 0.0, := 0.0;
var p188 >= 0.0, := 0.0;
var p189 >= 0.0, := 0.0;
var p190 >= 0.0, := 0.0;
var p191 >= 0.0, := 0.0;
var p192 >= 0.0, := 0.0;
var p193 >= 0.0, := 0.0;
var p194 >= 0.0, := 0.0;
var p195 >= 0.0, := 0.0;
var p196 >= 0.0, := 0.0;
var p197 >= 0.0, := 0.0;
var p198 >= 0.0, := 0.0;
var p199 >= 0.0, := 0.0;
var p200 >= 0.0, := 0.0;
var p201 >= 0.0, := 0.0;
var p202 >= 0.0, := 0.0;
var p203 >= 0.0, := 0.0;
var p204 >= 0.0, := 0.0;
var p205 >= 0.0, := 0.0;
var p206 >= 0.0, := 0.0;
var p207 >= 0.0, := 0.0;
var p208 >= 0.0, := 0.0;
var p209 >= 0.0, := 0.0;
var p210 >= 0.0, := 0.0;
var p211 >= 0.0, := 0.0;
var p212 >= 0.0, := 0.0;
var p213 >= 0.0, := 0.0;
var p214 >= 0.0, := 0.0;
var p215 >= 0.0, := 0.0;
var p216 >= 0.0, := 0.0;
var p217 >= 0.0, := 0.0;
var p218 >= 0.0, := 0.0;
var p219 >= 0.0, := 0.0;
var p220 >= 0.0, := 0.0;
var p221 >= 0.0, := 0.0;
var p222 >= 0.0, := 0.0;
var p223 >= 0.0, := 0.0;
var p224 >= 0.0, := 0.0;
var p225 >= 0.0, := 0.0;
var p226 >= 0.0, := 0.0;
var p227 >= 0.0, := 0.0;
var p228 >= 0.0, := 0.0;
var p229 >= 0.0, := 0.0;
var p230 >= 0.0, := 0.0;
var p231 >= 0.0, := 0.0;
var p232 >= 0.0, := 0.0;
var p233 >= 0.0, := 0.0;
var p234 >= 0.0, := 0.0;
var p235 >= 0.0, := 0.0;
var p236 >= 0.0, := 0.0;
var p237 >= 0.0, := 0.0;
var p238 >= 0.0, := 0.0;
var p239 >= 0.0, := 0.0;
var p240 >= 0.0, := 0.0;
var p241 >= 0.0, := 0.0;
var p242 >= 0.0, := 0.0;
var p243 >= 0.0, := 0.0;
var p244 >= 0.0, := 0.0;
var p245 >= 0.0, := 0.0;
var p246 >= 0.0, := 0.0;
var p247 >= 0.0, := 0.0;
var p248 >= 0.0, := 0.0;
var p249 >= 0.0, := 0.0;
var p250 >= 0.0, := 0.0;
var p251 >= 0.0, := 0.0;
var p252 >= 0.0, := 0.0;
var p253 >= 0.0, := 0.0;
var p254 >= 0.0, := 0.0;
var p255 >= 0.0, := 0.0;
var a1_1 >= 0.0 , := 0.0 ,  <= 1.0;
var a1_2 >= 0.0 , := 0.0 ,  <= 1.0;
var a1_3 >= 0.0 , := 0.0 ,  <= 1.0;
var a1_4 >= 0.0 , := 0.0 ,  <= 1.0;
var a1_5 >= 0.0 , := 0.0 ,  <= 1.0;
var a1_6 >= 0.0 , := 0.0 ,  <= 1.0;
var a1_7 >= 0.0 , := 0.0 ,  <= 1.0;
var a1_8 >= 0.0 , := 0.0 ,  <= 1.0;
var a2_1 >= 0.0 , := 0.0 ,  <= 1.0;
var a2_2 >= 0.0 , := 0.0 ,  <= 1.0;
var a2_3 >= 0.0 , := 0.0 ,  <= 1.0;
var a2_4 >= 0.0 , := 0.0 ,  <= 1.0;
var a2_5 >= 0.0 , := 0.0 ,  <= 1.0;
var a2_6 >= 0.0 , := 0.0 ,  <= 1.0;
var a2_7 >= 0.0 , := 0.0 ,  <= 1.0;
var a2_8 >= 0.0 , := 0.0 ,  <= 1.0;
var a3_1 >= 0.0 , := 0.0 ,  <= 1.0;
var a3_2 >= 0.0 , := 0.0 ,  <= 1.0;
var a3_3 >= 0.0 , := 0.0 ,  <= 1.0;
var a3_4 >= 0.0 , := 0.0 ,  <= 1.0;
var a3_5 >= 0.0 , := 0.0 ,  <= 1.0;
var a3_6 >= 0.0 , := 0.0 ,  <= 1.0;
var a3_7 >= 0.0 , := 0.0 ,  <= 1.0;
var a3_8 >= 0.0 , := 0.0 ,  <= 1.0;
var a4_1 >= 0.0 , := 0.0 ,  <= 1.0;
var a4_2 >= 0.0 , := 0.0 ,  <= 1.0;
var a4_3 >= 0.0 , := 0.0 ,  <= 1.0;
var a4_4 >= 0.0 , := 0.0 ,  <= 1.0;
var a4_5 >= 0.0 , := 0.0 ,  <= 1.0;
var a4_6 >= 0.0 , := 0.0 ,  <= 1.0;
var a4_7 >= 0.0 , := 0.0 ,  <= 1.0;
var a4_8 >= 0.0 , := 0.0 ,  <= 1.0;
var a5_1 >= 0.0 , := 0.0 ,  <= 1.0;
var a5_2 >= 0.0 , := 0.0 ,  <= 1.0;
var a5_3 >= 0.0 , := 0.0 ,  <= 1.0;
var a5_4 >= 0.0 , := 0.0 ,  <= 1.0;
var a5_5 >= 0.0 , := 0.0 ,  <= 1.0;
var a5_6 >= 0.0 , := 0.0 ,  <= 1.0;
var a5_7 >= 0.0 , := 0.0 ,  <= 1.0;
var a5_8 >= 0.0 , := 0.0 ,  <= 1.0;
var a6_1 >= 0.0 , := 0.0 ,  <= 1.0;
var a6_2 >= 0.0 , := 0.0 ,  <= 1.0;
var a6_3 >= 0.0 , := 0.0 ,  <= 1.0;
var a6_4 >= 0.0 , := 0.0 ,  <= 1.0;
var a6_5 >= 0.0 , := 0.0 ,  <= 1.0;
var a6_6 >= 0.0 , := 0.0 ,  <= 1.0;
var a6_7 >= 0.0 , := 0.0 ,  <= 1.0;
var a6_8 >= 0.0 , := 0.0 ,  <= 1.0;
var a7_1 >= 0.0 , := 0.0 ,  <= 1.0;
var a7_2 >= 0.0 , := 0.0 ,  <= 1.0;
var a7_3 >= 0.0 , := 0.0 ,  <= 1.0;
var a7_4 >= 0.0 , := 0.0 ,  <= 1.0;
var a7_5 >= 0.0 , := 0.0 ,  <= 1.0;
var a7_6 >= 0.0 , := 0.0 ,  <= 1.0;
var a7_7 >= 0.0 , := 0.0 ,  <= 1.0;
var a7_8 >= 0.0 , := 0.0 ,  <= 1.0;
var a8_1 >= 0.0 , := 0.0 ,  <= 1.0;
var a8_2 >= 0.0 , := 0.0 ,  <= 1.0;
var a8_3 >= 0.0 , := 0.0 ,  <= 1.0;
var a8_4 >= 0.0 , := 0.0 ,  <= 1.0;
var a8_5 >= 0.0 , := 0.0 ,  <= 1.0;
var a8_6 >= 0.0 , := 0.0 ,  <= 1.0;
var a8_7 >= 0.0 , := 0.0 ,  <= 1.0;
var a8_8 >= 0.0 , := 0.0 ,  <= 1.0;
var r1 >= 1.0d-6, := 0.1;
var c1 >= 1.0d-6, := 0.1;
var r2 >= 1.0d-6, := 0.1;
var c2 >= 1.0d-6, := 0.1;
var r3 >= 1.0d-6, := 0.1;
var c3 >= 1.0d-6, := 0.1;
var r4 >= 1.0d-6, := 0.1;
var c4 >= 1.0d-6, := 0.1;
var r5 >= 1.0d-6, := 0.1;
var c5 >= 1.0d-6, := 0.1;
var r6 >= 1.0d-6, := 0.1;
var c6 >= 1.0d-6, := 0.1;
var r7 >= 1.0d-6, := 0.1;
var c7 >= 1.0d-6, := 0.1;
var r8 >= 1.0d-6, := 0.1;
var c8 >= 1.0d-6, := 0.1;

minimize obj:
	(-exp((log ( r1 )  + log ( r2 )  + log ( r3 )  + log ( r4 )  + log ( r5 )  + log ( r6 )  + log ( r7 )  + log ( r8 ) ))) + (-exp((log ( c1 )  + log ( c2 )  + log ( c3 )  + log ( c4 )  + log ( c5 )  + log ( c6 )  + log ( c7 )  + log ( c8 ) ))) + p255;

subject to pe3:
	a2_7 * a1_8 + a2_8 * a1_7 - p3 = 0;
subject to pe5:
	a2_6 * a1_8 + a2_8 * a1_6 - p5 = 0;
subject to pe6:
	a2_6 * a1_7 + a2_7 * a1_6 - p6 = 0;
subject to pe7:
	a3_6 * p3 + a3_7 * p5 + a3_8 * p6 - p7 = 0;
subject to pe9:
	a2_5 * a1_8 + a2_8 * a1_5 - p9 = 0;
subject to pe10:
	a2_5 * a1_7 + a2_7 * a1_5 - p10 = 0;
subject to pe11:
	a3_5 * p3 + a3_7 * p9 + a3_8 * p10 - p11 = 0;
subject to pe12:
	a2_5 * a1_6 + a2_6 * a1_5 - p12 = 0;
subject to pe13:
	a3_5 * p5 + a3_6 * p9 + a3_8 * p12 - p13 = 0;
subject to pe14:
	a3_5 * p6 + a3_6 * p10 + a3_7 * p12 - p14 = 0;
subject to pe15:
	a4_5 * p7 + a4_6 * p11 + a4_7 * p13 + a4_8 * p14 - p15 = 0;
subject to pe17:
	a2_4 * a1_8 + a2_8 * a1_4 - p17 = 0;
subject to pe18:
	a2_4 * a1_7 + a2_7 * a1_4 - p18 = 0;
subject to pe19:
	a3_4 * p3 + a3_7 * p17 + a3_8 * p18 - p19 = 0;
subject to pe20:
	a2_4 * a1_6 + a2_6 * a1_4 - p20 = 0;
subject to pe21:
	a3_4 * p5 + a3_6 * p17 + a3_8 * p20 - p21 = 0;
subject to pe22:
	a3_4 * p6 + a3_6 * p18 + a3_7 * p20 - p22 = 0;
subject to pe23:
	a4_4 * p7 + a4_6 * p19 + a4_7 * p21 + a4_8 * p22 - p23 = 0;
subject to pe24:
	a2_4 * a1_5 + a2_5 * a1_4 - p24 = 0;
subject to pe25:
	a3_4 * p9 + a3_5 * p17 + a3_8 * p24 - p25 = 0;
subject to pe26:
	a3_4 * p10 + a3_5 * p18 + a3_7 * p24 - p26 = 0;
subject to pe27:
	a4_4 * p11 + a4_5 * p19 + a4_7 * p25 + a4_8 * p26 - p27 = 0;
subject to pe28:
	a3_4 * p12 + a3_5 * p20 + a3_6 * p24 - p28 = 0;
subject to pe29:
	a4_4 * p13 + a4_5 * p21 + a4_6 * p25 + a4_8 * p28 - p29 = 0;
subject to pe30:
	a4_4 * p14 + a4_5 * p22 + a4_6 * p26 + a4_7 * p28 - p30 = 0;
subject to pe31:
	a5_4 * p15 + a5_5 * p23 + a5_6 * p27 + a5_7 * p29 + a5_8 * p30 - p31 = 0;
subject to pe33:
	a2_3 * a1_8 + a2_8 * a1_3 - p33 = 0;
subject to pe34:
	a2_3 * a1_7 + a2_7 * a1_3 - p34 = 0;
subject to pe35:
	a3_3 * p3 + a3_7 * p33 + a3_8 * p34 - p35 = 0;
subject to pe36:
	a2_3 * a1_6 + a2_6 * a1_3 - p36 = 0;
subject to pe37:
	a3_3 * p5 + a3_6 * p33 + a3_8 * p36 - p37 = 0;
subject to pe38:
	a3_3 * p6 + a3_6 * p34 + a3_7 * p36 - p38 = 0;
subject to pe39:
	a4_3 * p7 + a4_6 * p35 + a4_7 * p37 + a4_8 * p38 - p39 = 0;
subject to pe40:
	a2_3 * a1_5 + a2_5 * a1_3 - p40 = 0;
subject to pe41:
	a3_3 * p9 + a3_5 * p33 + a3_8 * p40 - p41 = 0;
subject to pe42:
	a3_3 * p10 + a3_5 * p34 + a3_7 * p40 - p42 = 0;
subject to pe43:
	a4_3 * p11 + a4_5 * p35 + a4_7 * p41 + a4_8 * p42 - p43 = 0;
subject to pe44:
	a3_3 * p12 + a3_5 * p36 + a3_6 * p40 - p44 = 0;
subject to pe45:
	a4_3 * p13 + a4_5 * p37 + a4_6 * p41 + a4_8 * p44 - p45 = 0;
subject to pe46:
	a4_3 * p14 + a4_5 * p38 + a4_6 * p42 + a4_7 * p44 - p46 = 0;
subject to pe47:
	a5_3 * p15 + a5_5 * p39 + a5_6 * p43 + a5_7 * p45 + a5_8 * p46 - p47 = 0;
subject to pe48:
	a2_3 * a1_4 + a2_4 * a1_3 - p48 = 0;
subject to pe49:
	a3_3 * p17 + a3_4 * p33 + a3_8 * p48 - p49 = 0;
subject to pe50:
	a3_3 * p18 + a3_4 * p34 + a3_7 * p48 - p50 = 0;
subject to pe51:
	a4_3 * p19 + a4_4 * p35 + a4_7 * p49 + a4_8 * p50 - p51 = 0;
subject to pe52:
	a3_3 * p20 + a3_4 * p36 + a3_6 * p48 - p52 = 0;
subject to pe53:
	a4_3 * p21 + a4_4 * p37 + a4_6 * p49 + a4_8 * p52 - p53 = 0;
subject to pe54:
	a4_3 * p22 + a4_4 * p38 + a4_6 * p50 + a4_7 * p52 - p54 = 0;
subject to pe55:
	a5_3 * p23 + a5_4 * p39 + a5_6 * p51 + a5_7 * p53 + a5_8 * p54 - p55 = 0;
subject to pe56:
	a3_3 * p24 + a3_4 * p40 + a3_5 * p48 - p56 = 0;
subject to pe57:
	a4_3 * p25 + a4_4 * p41 + a4_5 * p49 + a4_8 * p56 - p57 = 0;
subject to pe58:
	a4_3 * p26 + a4_4 * p42 + a4_5 * p50 + a4_7 * p56 - p58 = 0;
subject to pe59:
	a5_3 * p27 + a5_4 * p43 + a5_5 * p51 + a5_7 * p57 + a5_8 * p58 - p59 = 0;
subject to pe60:
	a4_3 * p28 + a4_4 * p44 + a4_5 * p52 + a4_6 * p56 - p60 = 0;
subject to pe61:
	a5_3 * p29 + a5_4 * p45 + a5_5 * p53 + a5_6 * p57 + a5_8 * p60 - p61 = 0;
subject to pe62:
	a5_3 * p30 + a5_4 * p46 + a5_5 * p54 + a5_6 * p58 + a5_7 * p60 - p62 = 0;
subject to pe63:
	a6_3 * p31 + a6_4 * p47 + a6_5 * p55 + a6_6 * p59 + a6_7 * p61 + a6_8 * p62 - p63 = 0;
subject to pe65:
	a2_2 * a1_8 + a2_8 * a1_2 - p65 = 0;
subject to pe66:
	a2_2 * a1_7 + a2_7 * a1_2 - p66 = 0;
subject to pe67:
	a3_2 * p3 + a3_7 * p65 + a3_8 * p66 - p67 = 0;
subject to pe68:
	a2_2 * a1_6 + a2_6 * a1_2 - p68 = 0;
subject to pe69:
	a3_2 * p5 + a3_6 * p65 + a3_8 * p68 - p69 = 0;
subject to pe70:
	a3_2 * p6 + a3_6 * p66 + a3_7 * p68 - p70 = 0;
subject to pe71:
	a4_2 * p7 + a4_6 * p67 + a4_7 * p69 + a4_8 * p70 - p71 = 0;
subject to pe72:
	a2_2 * a1_5 + a2_5 * a1_2 - p72 = 0;
subject to pe73:
	a3_2 * p9 + a3_5 * p65 + a3_8 * p72 - p73 = 0;
subject to pe74:
	a3_2 * p10 + a3_5 * p66 + a3_7 * p72 - p74 = 0;
subject to pe75:
	a4_2 * p11 + a4_5 * p67 + a4_7 * p73 + a4_8 * p74 - p75 = 0;
subject to pe76:
	a3_2 * p12 + a3_5 * p68 + a3_6 * p72 - p76 = 0;
subject to pe77:
	a4_2 * p13 + a4_5 * p69 + a4_6 * p73 + a4_8 * p76 - p77 = 0;
subject to pe78:
	a4_2 * p14 + a4_5 * p70 + a4_6 * p74 + a4_7 * p76 - p78 = 0;
subject to pe79:
	a5_2 * p15 + a5_5 * p71 + a5_6 * p75 + a5_7 * p77 + a5_8 * p78 - p79 = 0;
subject to pe80:
	a2_2 * a1_4 + a2_4 * a1_2 - p80 = 0;
subject to pe81:
	a3_2 * p17 + a3_4 * p65 + a3_8 * p80 - p81 = 0;
subject to pe82:
	a3_2 * p18 + a3_4 * p66 + a3_7 * p80 - p82 = 0;
subject to pe83:
	a4_2 * p19 + a4_4 * p67 + a4_7 * p81 + a4_8 * p82 - p83 = 0;
subject to pe84:
	a3_2 * p20 + a3_4 * p68 + a3_6 * p80 - p84 = 0;
subject to pe85:
	a4_2 * p21 + a4_4 * p69 + a4_6 * p81 + a4_8 * p84 - p85 = 0;
subject to pe86:
	a4_2 * p22 + a4_4 * p70 + a4_6 * p82 + a4_7 * p84 - p86 = 0;
subject to pe87:
	a5_2 * p23 + a5_4 * p71 + a5_6 * p83 + a5_7 * p85 + a5_8 * p86 - p87 = 0;
subject to pe88:
	a3_2 * p24 + a3_4 * p72 + a3_5 * p80 - p88 = 0;
subject to pe89:
	a4_2 * p25 + a4_4 * p73 + a4_5 * p81 + a4_8 * p88 - p89 = 0;
subject to pe90:
	a4_2 * p26 + a4_4 * p74 + a4_5 * p82 + a4_7 * p88 - p90 = 0;
subject to pe91:
	a5_2 * p27 + a5_4 * p75 + a5_5 * p83 + a5_7 * p89 + a5_8 * p90 - p91 = 0;
subject to pe92:
	a4_2 * p28 + a4_4 * p76 + a4_5 * p84 + a4_6 * p88 - p92 = 0;
subject to pe93:
	a5_2 * p29 + a5_4 * p77 + a5_5 * p85 + a5_6 * p89 + a5_8 * p92 - p93 = 0;
subject to pe94:
	a5_2 * p30 + a5_4 * p78 + a5_5 * p86 + a5_6 * p90 + a5_7 * p92 - p94 = 0;
subject to pe95:
	a6_2 * p31 + a6_4 * p79 + a6_5 * p87 + a6_6 * p91 + a6_7 * p93 + a6_8 * p94 - p95 = 0;
subject to pe96:
	a2_2 * a1_3 + a2_3 * a1_2 - p96 = 0;
subject to pe97:
	a3_2 * p33 + a3_3 * p65 + a3_8 * p96 - p97 = 0;
subject to pe98:
	a3_2 * p34 + a3_3 * p66 + a3_7 * p96 - p98 = 0;
subject to pe99:
	a4_2 * p35 + a4_3 * p67 + a4_7 * p97 + a4_8 * p98 - p99 = 0;
subject to pe100:
	a3_2 * p36 + a3_3 * p68 + a3_6 * p96 - p100 = 0;
subject to pe101:
	a4_2 * p37 + a4_3 * p69 + a4_6 * p97 + a4_8 * p100 - p101 = 0;
subject to pe102:
	a4_2 * p38 + a4_3 * p70 + a4_6 * p98 + a4_7 * p100 - p102 = 0;
subject to pe103:
	a5_2 * p39 + a5_3 * p71 + a5_6 * p99 + a5_7 * p101 + a5_8 * p102 - p103 = 0;
subject to pe104:
	a3_2 * p40 + a3_3 * p72 + a3_5 * p96 - p104 = 0;
subject to pe105:
	a4_2 * p41 + a4_3 * p73 + a4_5 * p97 + a4_8 * p104 - p105 = 0;
subject to pe106:
	a4_2 * p42 + a4_3 * p74 + a4_5 * p98 + a4_7 * p104 - p106 = 0;
subject to pe107:
	a5_2 * p43 + a5_3 * p75 + a5_5 * p99 + a5_7 * p105 + a5_8 * p106 - p107 = 0;
subject to pe108:
	a4_2 * p44 + a4_3 * p76 + a4_5 * p100 + a4_6 * p104 - p108 = 0;
subject to pe109:
	a5_2 * p45 + a5_3 * p77 + a5_5 * p101 + a5_6 * p105 + a5_8 * p108 - p109 = 0;
subject to pe110:
	a5_2 * p46 + a5_3 * p78 + a5_5 * p102 + a5_6 * p106 + a5_7 * p108 - p110 = 0;
subject to pe111:
	a6_2 * p47 + a6_3 * p79 + a6_5 * p103 + a6_6 * p107 + a6_7 * p109 + a6_8 * p110 - p111 = 0;
subject to pe112:
	a3_2 * p48 + a3_3 * p80 + a3_4 * p96 - p112 = 0;
subject to pe113:
	a4_2 * p49 + a4_3 * p81 + a4_4 * p97 + a4_8 * p112 - p113 = 0;
subject to pe114:
	a4_2 * p50 + a4_3 * p82 + a4_4 * p98 + a4_7 * p112 - p114 = 0;
subject to pe115:
	a5_2 * p51 + a5_3 * p83 + a5_4 * p99 + a5_7 * p113 + a5_8 * p114 - p115 = 0;
subject to pe116:
	a4_2 * p52 + a4_3 * p84 + a4_4 * p100 + a4_6 * p112 - p116 = 0;
subject to pe117:
	a5_2 * p53 + a5_3 * p85 + a5_4 * p101 + a5_6 * p113 + a5_8 * p116 - p117 = 0;
subject to pe118:
	a5_2 * p54 + a5_3 * p86 + a5_4 * p102 + a5_6 * p114 + a5_7 * p116 - p118 = 0;
subject to pe119:
	a6_2 * p55 + a6_3 * p87 + a6_4 * p103 + a6_6 * p115 + a6_7 * p117 + a6_8 * p118 - p119 = 0;
subject to pe120:
	a4_2 * p56 + a4_3 * p88 + a4_4 * p104 + a4_5 * p112 - p120 = 0;
subject to pe121:
	a5_2 * p57 + a5_3 * p89 + a5_4 * p105 + a5_5 * p113 + a5_8 * p120 - p121 = 0;
subject to pe122:
	a5_2 * p58 + a5_3 * p90 + a5_4 * p106 + a5_5 * p114 + a5_7 * p120 - p122 = 0;
subject to pe123:
	a6_2 * p59 + a6_3 * p91 + a6_4 * p107 + a6_5 * p115 + a6_7 * p121 + a6_8 * p122 - p123 = 0;
subject to pe124:
	a5_2 * p60 + a5_3 * p92 + a5_4 * p108 + a5_5 * p116 + a5_6 * p120 - p124 = 0;
subject to pe125:
	a6_2 * p61 + a6_3 * p93 + a6_4 * p109 + a6_5 * p117 + a6_6 * p121 + a6_8 * p124 - p125 = 0;
subject to pe126:
	a6_2 * p62 + a6_3 * p94 + a6_4 * p110 + a6_5 * p118 + a6_6 * p122 + a6_7 * p124 - p126 = 0;
subject to pe127:
	a7_2 * p63 + a7_3 * p95 + a7_4 * p111 + a7_5 * p119 + a7_6 * p123 + a7_7 * p125 + a7_8 * p126 - p127 = 0;
subject to pe129:
	a2_1 * a1_8 + a2_8 * a1_1 - p129 = 0;
subject to pe130:
	a2_1 * a1_7 + a2_7 * a1_1 - p130 = 0;
subject to pe131:
	a3_1 * p3 + a3_7 * p129 + a3_8 * p130 - p131 = 0;
subject to pe132:
	a2_1 * a1_6 + a2_6 * a1_1 - p132 = 0;
subject to pe133:
	a3_1 * p5 + a3_6 * p129 + a3_8 * p132 - p133 = 0;
subject to pe134:
	a3_1 * p6 + a3_6 * p130 + a3_7 * p132 - p134 = 0;
subject to pe135:
	a4_1 * p7 + a4_6 * p131 + a4_7 * p133 + a4_8 * p134 - p135 = 0;
subject to pe136:
	a2_1 * a1_5 + a2_5 * a1_1 - p136 = 0;
subject to pe137:
	a3_1 * p9 + a3_5 * p129 + a3_8 * p136 - p137 = 0;
subject to pe138:
	a3_1 * p10 + a3_5 * p130 + a3_7 * p136 - p138 = 0;
subject to pe139:
	a4_1 * p11 + a4_5 * p131 + a4_7 * p137 + a4_8 * p138 - p139 = 0;
subject to pe140:
	a3_1 * p12 + a3_5 * p132 + a3_6 * p136 - p140 = 0;
subject to pe141:
	a4_1 * p13 + a4_5 * p133 + a4_6 * p137 + a4_8 * p140 - p141 = 0;
subject to pe142:
	a4_1 * p14 + a4_5 * p134 + a4_6 * p138 + a4_7 * p140 - p142 = 0;
subject to pe143:
	a5_1 * p15 + a5_5 * p135 + a5_6 * p139 + a5_7 * p141 + a5_8 * p142 - p143 = 0;
subject to pe144:
	a2_1 * a1_4 + a2_4 * a1_1 - p144 = 0;
subject to pe145:
	a3_1 * p17 + a3_4 * p129 + a3_8 * p144 - p145 = 0;
subject to pe146:
	a3_1 * p18 + a3_4 * p130 + a3_7 * p144 - p146 = 0;
subject to pe147:
	a4_1 * p19 + a4_4 * p131 + a4_7 * p145 + a4_8 * p146 - p147 = 0;
subject to pe148:
	a3_1 * p20 + a3_4 * p132 + a3_6 * p144 - p148 = 0;
subject to pe149:
	a4_1 * p21 + a4_4 * p133 + a4_6 * p145 + a4_8 * p148 - p149 = 0;
subject to pe150:
	a4_1 * p22 + a4_4 * p134 + a4_6 * p146 + a4_7 * p148 - p150 = 0;
subject to pe151:
	a5_1 * p23 + a5_4 * p135 + a5_6 * p147 + a5_7 * p149 + a5_8 * p150 - p151 = 0;
subject to pe152:
	a3_1 * p24 + a3_4 * p136 + a3_5 * p144 - p152 = 0;
subject to pe153:
	a4_1 * p25 + a4_4 * p137 + a4_5 * p145 + a4_8 * p152 - p153 = 0;
subject to pe154:
	a4_1 * p26 + a4_4 * p138 + a4_5 * p146 + a4_7 * p152 - p154 = 0;
subject to pe155:
	a5_1 * p27 + a5_4 * p139 + a5_5 * p147 + a5_7 * p153 + a5_8 * p154 - p155 = 0;
subject to pe156:
	a4_1 * p28 + a4_4 * p140 + a4_5 * p148 + a4_6 * p152 - p156 = 0;
subject to pe157:
	a5_1 * p29 + a5_4 * p141 + a5_5 * p149 + a5_6 * p153 + a5_8 * p156 - p157 = 0;
subject to pe158:
	a5_1 * p30 + a5_4 * p142 + a5_5 * p150 + a5_6 * p154 + a5_7 * p156 - p158 = 0;
subject to pe159:
	a6_1 * p31 + a6_4 * p143 + a6_5 * p151 + a6_6 * p155 + a6_7 * p157 + a6_8 * p158 - p159 = 0;
subject to pe160:
	a2_1 * a1_3 + a2_3 * a1_1 - p160 = 0;
subject to pe161:
	a3_1 * p33 + a3_3 * p129 + a3_8 * p160 - p161 = 0;
subject to pe162:
	a3_1 * p34 + a3_3 * p130 + a3_7 * p160 - p162 = 0;
subject to pe163:
	a4_1 * p35 + a4_3 * p131 + a4_7 * p161 + a4_8 * p162 - p163 = 0;
subject to pe164:
	a3_1 * p36 + a3_3 * p132 + a3_6 * p160 - p164 = 0;
subject to pe165:
	a4_1 * p37 + a4_3 * p133 + a4_6 * p161 + a4_8 * p164 - p165 = 0;
subject to pe166:
	a4_1 * p38 + a4_3 * p134 + a4_6 * p162 + a4_7 * p164 - p166 = 0;
subject to pe167:
	a5_1 * p39 + a5_3 * p135 + a5_6 * p163 + a5_7 * p165 + a5_8 * p166 - p167 = 0;
subject to pe168:
	a3_1 * p40 + a3_3 * p136 + a3_5 * p160 - p168 = 0;
subject to pe169:
	a4_1 * p41 + a4_3 * p137 + a4_5 * p161 + a4_8 * p168 - p169 = 0;
subject to pe170:
	a4_1 * p42 + a4_3 * p138 + a4_5 * p162 + a4_7 * p168 - p170 = 0;
subject to pe171:
	a5_1 * p43 + a5_3 * p139 + a5_5 * p163 + a5_7 * p169 + a5_8 * p170 - p171 = 0;
subject to pe172:
	a4_1 * p44 + a4_3 * p140 + a4_5 * p164 + a4_6 * p168 - p172 = 0;
subject to pe173:
	a5_1 * p45 + a5_3 * p141 + a5_5 * p165 + a5_6 * p169 + a5_8 * p172 - p173 = 0;
subject to pe174:
	a5_1 * p46 + a5_3 * p142 + a5_5 * p166 + a5_6 * p170 + a5_7 * p172 - p174 = 0;
subject to pe175:
	a6_1 * p47 + a6_3 * p143 + a6_5 * p167 + a6_6 * p171 + a6_7 * p173 + a6_8 * p174 - p175 = 0;
subject to pe176:
	a3_1 * p48 + a3_3 * p144 + a3_4 * p160 - p176 = 0;
subject to pe177:
	a4_1 * p49 + a4_3 * p145 + a4_4 * p161 + a4_8 * p176 - p177 = 0;
subject to pe178:
	a4_1 * p50 + a4_3 * p146 + a4_4 * p162 + a4_7 * p176 - p178 = 0;
subject to pe179:
	a5_1 * p51 + a5_3 * p147 + a5_4 * p163 + a5_7 * p177 + a5_8 * p178 - p179 = 0;
subject to pe180:
	a4_1 * p52 + a4_3 * p148 + a4_4 * p164 + a4_6 * p176 - p180 = 0;
subject to pe181:
	a5_1 * p53 + a5_3 * p149 + a5_4 * p165 + a5_6 * p177 + a5_8 * p180 - p181 = 0;
subject to pe182:
	a5_1 * p54 + a5_3 * p150 + a5_4 * p166 + a5_6 * p178 + a5_7 * p180 - p182 = 0;
subject to pe183:
	a6_1 * p55 + a6_3 * p151 + a6_4 * p167 + a6_6 * p179 + a6_7 * p181 + a6_8 * p182 - p183 = 0;
subject to pe184:
	a4_1 * p56 + a4_3 * p152 + a4_4 * p168 + a4_5 * p176 - p184 = 0;
subject to pe185:
	a5_1 * p57 + a5_3 * p153 + a5_4 * p169 + a5_5 * p177 + a5_8 * p184 - p185 = 0;
subject to pe186:
	a5_1 * p58 + a5_3 * p154 + a5_4 * p170 + a5_5 * p178 + a5_7 * p184 - p186 = 0;
subject to pe187:
	a6_1 * p59 + a6_3 * p155 + a6_4 * p171 + a6_5 * p179 + a6_7 * p185 + a6_8 * p186 - p187 = 0;
subject to pe188:
	a5_1 * p60 + a5_3 * p156 + a5_4 * p172 + a5_5 * p180 + a5_6 * p184 - p188 = 0;
subject to pe189:
	a6_1 * p61 + a6_3 * p157 + a6_4 * p173 + a6_5 * p181 + a6_6 * p185 + a6_8 * p188 - p189 = 0;
subject to pe190:
	a6_1 * p62 + a6_3 * p158 + a6_4 * p174 + a6_5 * p182 + a6_6 * p186 + a6_7 * p188 - p190 = 0;
subject to pe191:
	a7_1 * p63 + a7_3 * p159 + a7_4 * p175 + a7_5 * p183 + a7_6 * p187 + a7_7 * p189 + a7_8 * p190 - p191 = 0;
subject to pe192:
	a2_1 * a1_2 + a2_2 * a1_1 - p192 = 0;
subject to pe193:
	a3_1 * p65 + a3_2 * p129 + a3_8 * p192 - p193 = 0;
subject to pe194:
	a3_1 * p66 + a3_2 * p130 + a3_7 * p192 - p194 = 0;
subject to pe195:
	a4_1 * p67 + a4_2 * p131 + a4_7 * p193 + a4_8 * p194 - p195 = 0;
subject to pe196:
	a3_1 * p68 + a3_2 * p132 + a3_6 * p192 - p196 = 0;
subject to pe197:
	a4_1 * p69 + a4_2 * p133 + a4_6 * p193 + a4_8 * p196 - p197 = 0;
subject to pe198:
	a4_1 * p70 + a4_2 * p134 + a4_6 * p194 + a4_7 * p196 - p198 = 0;
subject to pe199:
	a5_1 * p71 + a5_2 * p135 + a5_6 * p195 + a5_7 * p197 + a5_8 * p198 - p199 = 0;
subject to pe200:
	a3_1 * p72 + a3_2 * p136 + a3_5 * p192 - p200 = 0;
subject to pe201:
	a4_1 * p73 + a4_2 * p137 + a4_5 * p193 + a4_8 * p200 - p201 = 0;
subject to pe202:
	a4_1 * p74 + a4_2 * p138 + a4_5 * p194 + a4_7 * p200 - p202 = 0;
subject to pe203:
	a5_1 * p75 + a5_2 * p139 + a5_5 * p195 + a5_7 * p201 + a5_8 * p202 - p203 = 0;
subject to pe204:
	a4_1 * p76 + a4_2 * p140 + a4_5 * p196 + a4_6 * p200 - p204 = 0;
subject to pe205:
	a5_1 * p77 + a5_2 * p141 + a5_5 * p197 + a5_6 * p201 + a5_8 * p204 - p205 = 0;
subject to pe206:
	a5_1 * p78 + a5_2 * p142 + a5_5 * p198 + a5_6 * p202 + a5_7 * p204 - p206 = 0;
subject to pe207:
	a6_1 * p79 + a6_2 * p143 + a6_5 * p199 + a6_6 * p203 + a6_7 * p205 + a6_8 * p206 - p207 = 0;
subject to pe208:
	a3_1 * p80 + a3_2 * p144 + a3_4 * p192 - p208 = 0;
subject to pe209:
	a4_1 * p81 + a4_2 * p145 + a4_4 * p193 + a4_8 * p208 - p209 = 0;
subject to pe210:
	a4_1 * p82 + a4_2 * p146 + a4_4 * p194 + a4_7 * p208 - p210 = 0;
subject to pe211:
	a5_1 * p83 + a5_2 * p147 + a5_4 * p195 + a5_7 * p209 + a5_8 * p210 - p211 = 0;
subject to pe212:
	a4_1 * p84 + a4_2 * p148 + a4_4 * p196 + a4_6 * p208 - p212 = 0;
subject to pe213:
	a5_1 * p85 + a5_2 * p149 + a5_4 * p197 + a5_6 * p209 + a5_8 * p212 - p213 = 0;
subject to pe214:
	a5_1 * p86 + a5_2 * p150 + a5_4 * p198 + a5_6 * p210 + a5_7 * p212 - p214 = 0;
subject to pe215:
	a6_1 * p87 + a6_2 * p151 + a6_4 * p199 + a6_6 * p211 + a6_7 * p213 + a6_8 * p214 - p215 = 0;
subject to pe216:
	a4_1 * p88 + a4_2 * p152 + a4_4 * p200 + a4_5 * p208 - p216 = 0;
subject to pe217:
	a5_1 * p89 + a5_2 * p153 + a5_4 * p201 + a5_5 * p209 + a5_8 * p216 - p217 = 0;
subject to pe218:
	a5_1 * p90 + a5_2 * p154 + a5_4 * p202 + a5_5 * p210 + a5_7 * p216 - p218 = 0;
subject to pe219:
	a6_1 * p91 + a6_2 * p155 + a6_4 * p203 + a6_5 * p211 + a6_7 * p217 + a6_8 * p218 - p219 = 0;
subject to pe220:
	a5_1 * p92 + a5_2 * p156 + a5_4 * p204 + a5_5 * p212 + a5_6 * p216 - p220 = 0;
subject to pe221:
	a6_1 * p93 + a6_2 * p157 + a6_4 * p205 + a6_5 * p213 + a6_6 * p217 + a6_8 * p220 - p221 = 0;
subject to pe222:
	a6_1 * p94 + a6_2 * p158 + a6_4 * p206 + a6_5 * p214 + a6_6 * p218 + a6_7 * p220 - p222 = 0;
subject to pe223:
	a7_1 * p95 + a7_2 * p159 + a7_4 * p207 + a7_5 * p215 + a7_6 * p219 + a7_7 * p221 + a7_8 * p222 - p223 = 0;
subject to pe224:
	a3_1 * p96 + a3_2 * p160 + a3_3 * p192 - p224 = 0;
subject to pe225:
	a4_1 * p97 + a4_2 * p161 + a4_3 * p193 + a4_8 * p224 - p225 = 0;
subject to pe226:
	a4_1 * p98 + a4_2 * p162 + a4_3 * p194 + a4_7 * p224 - p226 = 0;
subject to pe227:
	a5_1 * p99 + a5_2 * p163 + a5_3 * p195 + a5_7 * p225 + a5_8 * p226 - p227 = 0;
subject to pe228:
	a4_1 * p100 + a4_2 * p164 + a4_3 * p196 + a4_6 * p224 - p228 = 0;
subject to pe229:
	a5_1 * p101 + a5_2 * p165 + a5_3 * p197 + a5_6 * p225 + a5_8 * p228 - p229 = 0;
subject to pe230:
	a5_1 * p102 + a5_2 * p166 + a5_3 * p198 + a5_6 * p226 + a5_7 * p228 - p230 = 0;
subject to pe231:
	a6_1 * p103 + a6_2 * p167 + a6_3 * p199 + a6_6 * p227 + a6_7 * p229 + a6_8 * p230 - p231 = 0;
subject to pe232:
	a4_1 * p104 + a4_2 * p168 + a4_3 * p200 + a4_5 * p224 - p232 = 0;
subject to pe233:
	a5_1 * p105 + a5_2 * p169 + a5_3 * p201 + a5_5 * p225 + a5_8 * p232 - p233 = 0;
subject to pe234:
	a5_1 * p106 + a5_2 * p170 + a5_3 * p202 + a5_5 * p226 + a5_7 * p232 - p234 = 0;
subject to pe235:
	a6_1 * p107 + a6_2 * p171 + a6_3 * p203 + a6_5 * p227 + a6_7 * p233 + a6_8 * p234 - p235 = 0;
subject to pe236:
	a5_1 * p108 + a5_2 * p172 + a5_3 * p204 + a5_5 * p228 + a5_6 * p232 - p236 = 0;
subject to pe237:
	a6_1 * p109 + a6_2 * p173 + a6_3 * p205 + a6_5 * p229 + a6_6 * p233 + a6_8 * p236 - p237 = 0;
subject to pe238:
	a6_1 * p110 + a6_2 * p174 + a6_3 * p206 + a6_5 * p230 + a6_6 * p234 + a6_7 * p236 - p238 = 0;
subject to pe239:
	a7_1 * p111 + a7_2 * p175 + a7_3 * p207 + a7_5 * p231 + a7_6 * p235 + a7_7 * p237 + a7_8 * p238 - p239 = 0;
subject to pe240:
	a4_1 * p112 + a4_2 * p176 + a4_3 * p208 + a4_4 * p224 - p240 = 0;
subject to pe241:
	a5_1 * p113 + a5_2 * p177 + a5_3 * p209 + a5_4 * p225 + a5_8 * p240 - p241 = 0;
subject to pe242:
	a5_1 * p114 + a5_2 * p178 + a5_3 * p210 + a5_4 * p226 + a5_7 * p240 - p242 = 0;
subject to pe243:
	a6_1 * p115 + a6_2 * p179 + a6_3 * p211 + a6_4 * p227 + a6_7 * p241 + a6_8 * p242 - p243 = 0;
subject to pe244:
	a5_1 * p116 + a5_2 * p180 + a5_3 * p212 + a5_4 * p228 + a5_6 * p240 - p244 = 0;
subject to pe245:
	a6_1 * p117 + a6_2 * p181 + a6_3 * p213 + a6_4 * p229 + a6_6 * p241 + a6_8 * p244 - p245 = 0;
subject to pe246:
	a6_1 * p118 + a6_2 * p182 + a6_3 * p214 + a6_4 * p230 + a6_6 * p242 + a6_7 * p244 - p246 = 0;
subject to pe247:
	a7_1 * p119 + a7_2 * p183 + a7_3 * p215 + a7_4 * p231 + a7_6 * p243 + a7_7 * p245 + a7_8 * p246 - p247 = 0;
subject to pe248:
	a5_1 * p120 + a5_2 * p184 + a5_3 * p216 + a5_4 * p232 + a5_5 * p240 - p248 = 0;
subject to pe249:
	a6_1 * p121 + a6_2 * p185 + a6_3 * p217 + a6_4 * p233 + a6_5 * p241 + a6_8 * p248 - p249 = 0;
subject to pe250:
	a6_1 * p122 + a6_2 * p186 + a6_3 * p218 + a6_4 * p234 + a6_5 * p242 + a6_7 * p248 - p250 = 0;
subject to pe251:
	a7_1 * p123 + a7_2 * p187 + a7_3 * p219 + a7_4 * p235 + a7_5 * p243 + a7_7 * p249 + a7_8 * p250 - p251 = 0;
subject to pe252:
	a6_1 * p124 + a6_2 * p188 + a6_3 * p220 + a6_4 * p236 + a6_5 * p244 + a6_6 * p248 - p252 = 0;
subject to pe253:
	a7_1 * p125 + a7_2 * p189 + a7_3 * p221 + a7_4 * p237 + a7_5 * p245 + a7_6 * p249 + a7_8 * p252 - p253 = 0;
subject to pe254:
	a7_1 * p126 + a7_2 * p190 + a7_3 * p222 + a7_4 * p238 + a7_5 * p246 + a7_6 * p250 + a7_7 * p252 - p254 = 0;
subject to pe255:
	a8_1 * p127 + a8_2 * p191 + a8_3 * p223 + a8_4 * p239 + a8_5 * p247 + a8_6 * p251 + a8_7 * p253 + a8_8 * p254 - p255 = 0;
subject to r1_cons:
	a1_1 + a1_2 + a1_3 + a1_4 + a1_5 + a1_6 + a1_7 + a1_8 - r1 = 0;
subject to c1_cons:
	a1_1 + a2_1 + a3_1 + a4_1 + a5_1 + a6_1 + a7_1 + a8_1 - c1 = 0;
subject to c2_cons:
	a1_2 + a2_2 + a3_2 + a4_2 + a5_2 + a6_2 + a7_2 + a8_2 - c2 = 0;
subject to c3_cons:
	a1_3 + a2_3 + a3_3 + a4_3 + a5_3 + a6_3 + a7_3 + a8_3 - c3 = 0;
subject to c4_cons:
	a1_4 + a2_4 + a3_4 + a4_4 + a5_4 + a6_4 + a7_4 + a8_4 - c4 = 0;
subject to c5_cons:
	a1_5 + a2_5 + a3_5 + a4_5 + a5_5 + a6_5 + a7_5 + a8_5 - c5 = 0;
subject to c6_cons:
	a1_6 + a2_6 + a3_6 + a4_6 + a5_6 + a6_6 + a7_6 + a8_6 - c6 = 0;
subject to c7_cons:
	a1_7 + a2_7 + a3_7 + a4_7 + a5_7 + a6_7 + a7_7 + a8_7 - c7 = 0;
subject to c8_cons:
	a1_8 + a2_8 + a3_8 + a4_8 + a5_8 + a6_8 + a7_8 + a8_8 - c8 = 0;
subject to r2_cons:
	a2_1 + a2_2 + a2_3 + a2_4 + a2_5 + a2_6 + a2_7 + a2_8 - r2 = 0;
subject to r3_cons:
	a3_1 + a3_2 + a3_3 + a3_4 + a3_5 + a3_6 + a3_7 + a3_8 - r3 = 0;
subject to r4_cons:
	a4_1 + a4_2 + a4_3 + a4_4 + a4_5 + a4_6 + a4_7 + a4_8 - r4 = 0;
subject to r5_cons:
	a5_1 + a5_2 + a5_3 + a5_4 + a5_5 + a5_6 + a5_7 + a5_8 - r5 = 0;
subject to r6_cons:
	a6_1 + a6_2 + a6_3 + a6_4 + a6_5 + a6_6 + a6_7 + a6_8 - r6 = 0;
subject to r7_cons:
	a7_1 + a7_2 + a7_3 + a7_4 + a7_5 + a7_6 + a7_7 + a7_8 - r7 = 0;
subject to r8_cons:
	a8_1 + a8_2 + a8_3 + a8_4 + a8_5 + a8_6 + a8_7 + a8_8 - r8 = 0;
subject to sum_cons:
	r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8 - 8.0 = 0;

solve;
display p3;
display p5;
display p6;
display p7;
display p9;
display p10;
display p11;
display p12;
display p13;
display p14;
display p15;
display p17;
display p18;
display p19;
display p20;
display p21;
display p22;
display p23;
display p24;
display p25;
display p26;
display p27;
display p28;
display p29;
display p30;
display p31;
display p33;
display p34;
display p35;
display p36;
display p37;
display p38;
display p39;
display p40;
display p41;
display p42;
display p43;
display p44;
display p45;
display p46;
display p47;
display p48;
display p49;
display p50;
display p51;
display p52;
display p53;
display p54;
display p55;
display p56;
display p57;
display p58;
display p59;
display p60;
display p61;
display p62;
display p63;
display p65;
display p66;
display p67;
display p68;
display p69;
display p70;
display p71;
display p72;
display p73;
display p74;
display p75;
display p76;
display p77;
display p78;
display p79;
display p80;
display p81;
display p82;
display p83;
display p84;
display p85;
display p86;
display p87;
display p88;
display p89;
display p90;
display p91;
display p92;
display p93;
display p94;
display p95;
display p96;
display p97;
display p98;
display p99;
display p100;
display p101;
display p102;
display p103;
display p104;
display p105;
display p106;
display p107;
display p108;
display p109;
display p110;
display p111;
display p112;
display p113;
display p114;
display p115;
display p116;
display p117;
display p118;
display p119;
display p120;
display p121;
display p122;
display p123;
display p124;
display p125;
display p126;
display p127;
display p129;
display p130;
display p131;
display p132;
display p133;
display p134;
display p135;
display p136;
display p137;
display p138;
display p139;
display p140;
display p141;
display p142;
display p143;
display p144;
display p145;
display p146;
display p147;
display p148;
display p149;
display p150;
display p151;
display p152;
display p153;
display p154;
display p155;
display p156;
display p157;
display p158;
display p159;
display p160;
display p161;
display p162;
display p163;
display p164;
display p165;
display p166;
display p167;
display p168;
display p169;
display p170;
display p171;
display p172;
display p173;
display p174;
display p175;
display p176;
display p177;
display p178;
display p179;
display p180;
display p181;
display p182;
display p183;
display p184;
display p185;
display p186;
display p187;
display p188;
display p189;
display p190;
display p191;
display p192;
display p193;
display p194;
display p195;
display p196;
display p197;
display p198;
display p199;
display p200;
display p201;
display p202;
display p203;
display p204;
display p205;
display p206;
display p207;
display p208;
display p209;
display p210;
display p211;
display p212;
display p213;
display p214;
display p215;
display p216;
display p217;
display p218;
display p219;
display p220;
display p221;
display p222;
display p223;
display p224;
display p225;
display p226;
display p227;
display p228;
display p229;
display p230;
display p231;
display p232;
display p233;
display p234;
display p235;
display p236;
display p237;
display p238;
display p239;
display p240;
display p241;
display p242;
display p243;
display p244;
display p245;
display p246;
display p247;
display p248;
display p249;
display p250;
display p251;
display p252;
display p253;
display p254;
display p255;
display a1_1;
display a1_2;
display a1_3;
display a1_4;
display a1_5;
display a1_6;
display a1_7;
display a1_8;
display a2_1;
display a2_2;
display a2_3;
display a2_4;
display a2_5;
display a2_6;
display a2_7;
display a2_8;
display a3_1;
display a3_2;
display a3_3;
display a3_4;
display a3_5;
display a3_6;
display a3_7;
display a3_8;
display a4_1;
display a4_2;
display a4_3;
display a4_4;
display a4_5;
display a4_6;
display a4_7;
display a4_8;
display a5_1;
display a5_2;
display a5_3;
display a5_4;
display a5_5;
display a5_6;
display a5_7;
display a5_8;
display a6_1;
display a6_2;
display a6_3;
display a6_4;
display a6_5;
display a6_6;
display a6_7;
display a6_8;
display a7_1;
display a7_2;
display a7_3;
display a7_4;
display a7_5;
display a7_6;
display a7_7;
display a7_8;
display a8_1;
display a8_2;
display a8_3;
display a8_4;
display a8_5;
display a8_6;
display a8_7;
display a8_8;
display r1;
display c1;
display r2;
display c2;
display r3;
display c3;
display r4;
display c4;
display r5;
display c5;
display r6;
display c6;
display r7;
display c7;
display r8;
display c8;
display obj;
