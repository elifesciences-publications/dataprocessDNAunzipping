function [zs,workfjc,workwlc,ftot,zss,zds] = getStretchingWork(params)
% get the work required to to stretch a ss and a ds chain up to specific
% fractional extensions listed in zs
% and the appropriate fractional extensions of the ss and ds components of
% a hybrid chain held at forces given in ftot.

bkT = params.bss/params.kT; lp = params.lp; K = params.K; Lds = params.Lds; Lss = params.Lss; Kss = params.Kss; lpkT = lp/params.kT;

% get fractional extension of ss and ds chain for a range of forces
nf = 1000;
ftot = linspace(0.001,50,nf);
[zss,zds] = getChainExtFromForce(ftot, params);


%% get force-extension work for each total fractional extension
nz=1000; 
zs = linspace(1e-3,1.5,nz);
dz = zs(2)-zs(1);
fwlc = getWLCapproxF(zs,lpkT,K);
workwlc = cumsum(fwlc)*dz*2*Lds;

ffjc = zeros(1,nz);
for zc = 1:nz;
    ffjc(zc) = fzero(@(f) getFJCapproxZ(f,bkT,K) - zs(zc),10);
end
workfjc = cumsum(ffjc)*dz;