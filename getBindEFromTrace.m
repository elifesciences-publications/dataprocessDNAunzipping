function [fitbindE,zs,workfjc,workwlc,ftot,zss,zds] = getBindEFromTrace(binx,fsmooth, params)
% get binding energy between each DNA basepair from an averaged unzipping
% trace
% assumes ssDNA parameters have been previously fitted
% binx and fsmooth are smoothed distances and forces for the experimental
% curve
% params are obtained from ssdsForceExt.m
% returns energy in pN-nm (per bp)

[zs,workfjc,workwlc,ftot,zss,zds] = getStretchingWork(params);

%%
% for each extension, find the ss length that gives the smoothed force

% ss and ds fractional lengths for each force in fsmooth
zsx = interp1(ftot,zss,fsmooth);
zdx = interp1(ftot,zds,fsmooth);
Lsslist = (binx - zdx*2*params.Lds)./zsx;

%% work to unzip each bp

xind = find(binx>600 & Lsslist/params.ssdel/2<=params.Nss+1 );
%xind = xind:length(binx);
nunzip = Lsslist(xind)/params.ssdel/2;

% extensions for each integer n unzipped
nlist = 0:params.Nss;
extn = interp1(nunzip,binx(xind),nlist);
zssn = interp1(nunzip,zsx(xind),nlist);
zdsn = interp1(nunzip,zdx(xind),nlist);


nint = 100;
clear workdone fitbindE
for cc = 2:length(extn)    
    %%
    f0 = interp1(binx,fsmooth,extn(cc-1));
    wint = interp1(zs,workfjc,zssn(cc));
    workss2(cc-1) = 2*params.ssdel*wint;    
    
    fitbindE(cc-1) = f0*zssn(cc)*params.ssdel*2-workss2(cc-1);    
end

end