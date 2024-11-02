Uzip CalTopoEvol1.1.zip, then get two files: CalTopoEvol.rar and CalTopoEvol.wl 
Decompress CalTopoEvol.rar and get CalTopoEvol.mx
Open a Mathematica notebook and run the command: Get["Dir\CalTopoEvol.wl"], where Dir denotes the directory containing CalTopoEvol.wl
An example: 
CalTopoEvol["126.386",{{"\\[CapitalGamma][1,1]", {{{1}, 2, 11}, {{2}, 2, 8}, {{3}, 2, 10}, {{4}, 2, 5}}},
     {"M[1,1]", {{{1}, 4, 17}}},
     {"Z[1,1]", {{{1}, 2, 11}, {{2}, 2, 10}, {{3}, 2, 7}, {{4}, 2, 6}}},
     {"A[1,1]", {{{1}, 4, 17}}},
     {"R[1,1]", {{{1, 2}, 4, 17}}},
     {"X[1,1]", {{{1}, 2, 17}, {{2}, 2, 17}}}}, 1]
 where the first argument: "126.386" denotes name of the magnetic space group (in the BNS notation), the third argument so=1 denotes that SOC is 
 considered (when SOC can be neglected, so=0),
 the second argument denotes the high-symmetry point symmetry data.
 
 
