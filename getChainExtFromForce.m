function [zss,zds] = getChainExtFromForce(flist, options)
% for each of the forces in flist,
% calculate the fractional extension of single stranded and double stranded
% chain that gives you that force

% ds persistence length (in nm)
opt.lp = 35.4;
opt.K = 1020; % double stranded stretch modulus (in pN)
opt.Kss = 1000; %single stranded stretch modulus
opt.bss = 1.0348; %single stranded stiffness (in pN-nm);

opt.kT = 4.1;


% copy over input parameters
opt = copyStruct(options,opt,'addnew',1);

lpkT = opt.lp/opt.kT;
bkT = opt.bss/opt.kT;

zss = getFJCapproxZ(flist,bkT,opt.Kss);

% solve for force for a given WLC extension
for fc = 1:length(flist)    
    zds(fc) = fzero(@(x) getWLCapproxF(x,lpkT,opt.K)-flist(fc),0.9);
end


end