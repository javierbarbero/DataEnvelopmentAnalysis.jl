using COSMO, DataEnvelopmentAnalysis

X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];
Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];



test_normal = dea(X,Y, optimizer = DEAOptimizer(COSMO.Optimizer), names = dmus_names)
test_big = deabigdata(X, Y, optimizer = DEAOptimizer(COSMO.Optimizer))


