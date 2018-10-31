function [Evals, Estretch,forces,zplot,nlist] = getElandscape(bindE,params)
% get overall energy landscape for a given binding energy

% output:
% bindE = binding energies in kT per bp
% Evals = landscape of overall energies
% Estretch = landscape of stretching energies only (excludes binding
% energies)
% forces = landscape of forces corresponding to each z extension and number
% unzipped
% assumes parameters have been established

Lss = params.Lss; Lds = params.Lds;
nbp = params.Nss;

% zero out poorly defined unzipping energy values 
bindE(isnan(bindE)) = 0;
%  start = find(~isnan(bindE),1,'first');
%  if (~isempty(start) & start > 1)
%      bindE(1:start-1) = 0;
%  end

%% calculate work to stretch to fractional extensions
[zs,workfjc,workwlc,ftot,zss,zds] = getStretchingWork(params);
%%
% for each length unzipped, for each end distance
% calculate the corresponding ss and ds end distance and the work
nzplot = 200;
zplot = linspace(400,1100,nzplot);

nlist = [1:nbp, nbp+2];
sszval = zeros(nbp,nzplot);
dszval = sszval;
Wwlc = sszval;
Wfjc = sszval;
Evals = sszval; forces = sszval;
Estretch = Evals;
for cc = 1:length(nlist)
    nc = nlist(cc);
        Lss = nc*2*params.ssdel;
    ztot = zss*Lss+zds*2*Lds;
    sszval(cc,:)= interp1(ztot,zss, zplot);
    dszval(cc,:) = interp1(ztot,zds,zplot);
    
    forces(cc,:) = interp1(ztot,ftot,zplot);
    Wwlc(cc,:) = interp1(zs,workwlc,dszval(cc,:));
    Wfjc(cc,:) = interp1(zs,workfjc,sszval(cc,:));
    % energy in units of kT
    Evals(cc,:) = Wwlc(cc,:)/params.kT + Wfjc(cc,:)*(nc*2*params.ssdel)/params.kT + sum(bindE(1:min(nc,nbp)));  
    
    % stretching energy in units of kT only, excluding binding energy
    Estretch(cc,:) = Wwlc(cc,:)/params.kT + Wfjc(cc,:)*(nc*2*params.ssdel)/params.kT;
end

if (params.dodisplay>1)
    figure
    pcolor(zplot,nlist,Evals)
    shading flat
    xlabel('end to end distance (nm)')
    ylabel('number of bp unwrapped')
    colorbar
    title('energy in kT')
end
