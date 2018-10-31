function [seqshift,fitbindE,params] = getBindEfitSSParam(distvals,fvals,params)
% for a given set of force-extension data (interpolated to equispaced
% extensions), fit the appropriate ssDNA parameters
% 
% return the binding energy profile and the shift between corresponding
% features in the two halves of that profile
%
% INPUT:
% distinterp = equispaced extension distances (in nm), with an arbitrary
% unknown shift
% fvals = forces for each of the extensions
% binding energy profile obtained from two sequence copies
% params = ds and ss DNA params, excluding b and del for the ssDNA
%
% RETURNS:
% fitbindE = binding energies for all the basepairs
% params = parameter structure including bss (ssDNA Kuhn length) and ssdel
% (length in nm per bp of ssDNA)
% seqshift = shift between features in the two halves of the binding energy
% profile (in units of bp)

if (~isfield(params,'dodisplay'))
    params.dodisplay=0;        
end
if (~isfield(params,'fix'))
    params.fix=0;    
end

% fit after the last local minimum in the data
diffs = diff(fvals);
minind = find(diffs<0,1,'last')+1;
fitind = minind:length(distvals);

if (~params.fix)
    %% fit the ssDNA parameters    
    
    params.minf = fvals(fitind(1))-5;
    params.maxf = fvals(fitind(end))+40;
    
    % fit ss stiffness value
    
    cfit = nlinfit(distvals(fitind),fvals(fitind),@(c,x) ssdsForceExt(x,c(1)/params.kT,c(2),params),[1,params.Nsstot*0.57]);
    
    % ss stiffness in pN-nm
    params.bss = cfit(1);    
    params.ssdel = cfit(2)/params.Nsstot;
    params.Lss = cfit(2);
end

bkT = params.bss/params.kT;
%%
zext = linspace(distvals(fitind(1)),distvals(fitind(end)),50);
fmodel = ssdsForceExt(zext,bkT,params.Lss,params);

if (params.dodisplay>1)
    figure
    plot(distvals(fitind),fvals(fitind),'.',zext,fmodel,'LineWidth',2)
    xlabel('extension')
    ylabel('force')
    title('fit for fully unzipped region')
end

[fmodel] = ssdsForceExt(zext,params.bss/params.kT,params.Lss,params);


%% extract binding energy using these ssDNA parameters
nlist = 0:params.Nss;
[fitbindE,zs,workfjc,workwlc,ftot,zss,zds] = getBindEFromTrace(distvals,fvals,params);
fitbindE(isnan(fitbindE)) = 0;

% align two halves of sequence to see shift
ind1 = params.start1:params.end1;
ind2 = params.start2:params.end2;

s1 = fitbindE(ind1)-mean(fitbindE(ind1));
s2 = fitbindE(ind2)-mean(fitbindE(ind2));
tdel = finddelay(s2,s1);

% original shift between two pieces of the trace;
dn = nlist(ind1(1))-nlist(ind2(1));


if (params.dodisplay>1)
    % plot overlapping halves of the profile to check shift makes sense
    figure
    plot(nlist(ind1),fitbindE(ind1)/params.kT)
    hold all
    if (tdel<0)
        plot(nlist(ind2(1:end-abs(tdel)))+dn,fitbindE(ind2(abs(tdel)+1:end))/params.kT)
    else
        plot(nlist(ind2(abs(tdel)+1:end))+dn,fitbindE(ind2(1:end-abs(tdel)))/params.kT)
    end
    hold off
end

% overall shift between repeat sequences; should be 197
seqshift = abs(tdel+dn);
