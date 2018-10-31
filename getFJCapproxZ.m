function zvals = getFJCapproxZ(fvals,bkT,K)
% for an extensible freely joined chain, get the extension for different
% forces
% bkT = segment length /kT
% K = stretch modulus in force units
% zvals is given as a fraction of total length

zvals = (coth(fvals*bkT) - 1./(fvals*bkT)).*(1+fvals/K);
