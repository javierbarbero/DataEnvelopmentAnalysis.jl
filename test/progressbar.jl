using DataEnvelopmentAnalysis

# Cannot run the tests due to Ipot and GLPK not supported on ARM 64

X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];
Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];
optimizer = DEAOptimizer(:COSMO)
dea(X, Y, orient = :Input, rts = :CRS, optimizer = optimizer)