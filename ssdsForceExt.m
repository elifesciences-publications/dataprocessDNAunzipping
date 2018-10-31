function [fmodel] = ssdsForceExt(zext,bkT,Lss,options)
% get force for the given extensions
% of a chain composed of a ss region and 2 ds handles
% parameters:
% lp = ds persistence length
% kT = 4.1 pN nm
% K = ds stretch modulus
% Lds = ds length for one handle in nm
% Lss = total ss length in nm
% Kss = single stranded stretch modulus

% default values
opt.kT = 4.1;
opt.nf = 100;

% min and max force to use (in pN)
opt.minf = 0.1; 
opt.maxf = 30;

% copy over input parameters
opt = copyStruct(options,opt,'addnew',1);

lpkT = opt.lp/opt.kT;
Kss = opt.Kss;
Lds = opt.Lds;
K = opt.K;

ftot = linspace(opt.minf,opt.maxf,opt.nf);
zss = getFJCapproxZ(ftot,bkT,Kss);

% solve for force for a given WLC extension
for fc = 1:length(ftot)
    zds(fc) = fzero(@(x) getWLCapproxF(x,lpkT,K)-ftot(fc),0.9);
end

ztot = zss*Lss+zds*2*Lds;

% get forces corresponding to measured distances
fmodel = interp1(ztot,ftot,zext);

end
