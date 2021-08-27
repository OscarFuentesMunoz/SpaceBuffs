addpath(genpath(pwd))
format longG

OE1 = [1 0.1 0.1 1 1 1];
OE2 = [1.5 0.3 0.5 .5 .5 .5];

% Order of orbit elements in vector:
% sma - ecc - node - inclination - arg. per

MOID1 = MOID_SDG_win(OE1([1 2 4 3 5]), OE2([1 2 4 3 5]))


