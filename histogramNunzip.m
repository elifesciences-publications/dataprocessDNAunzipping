function [freq,binmid,cts] = histogramNunzip(zvals,fvals,nlist,params)
%% calculate normalized histogram for number of bp unzipped 
% from a given cluster of fixed-trap data points
% zvals = extension in nm
% fvals = force in pN
% nlist = list of # unzipped for which to calculate frequency values
% params = mechanical parameters for ss and ds DNA
% returns normalized frequency and bin positions
% also returns total (unnormalized) cts

nfs = 200;
fsample = linspace(min(fvals),max(fvals),nfs);
% fractional extensions for each force
[zss,zds] = getChainExtFromForce(fsample, params);


zsscur = interp1(fsample,zss,fvals);
zdscur = interp1(fsample,zds,fvals);

% single stranded length for each force
Lsslist = (zvals - zdscur*2*params.Lds)./zsscur;
% number of unzipped basepairs for each force
Nsslist = Lsslist/(params.ssdel*2);
%
binmid = nlist+0.5;
cts = hist(Nsslist,binmid);
normfact = (sum(cts) - 0.5*(cts(1)+cts(end)))*(binmid(2)-binmid(1));
freq = cts/normfact;
