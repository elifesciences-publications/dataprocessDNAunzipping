%% Run all processing steps to align extension based on an averaged
% bare DNA unzipping curve
% and then extract binding energies from hopping data
% use same ssDNA parameters obtained by fitting averaged bare DNA curves

dirname = '../data/hopping/Hopping/WT/';
hopname = '180228_067';
unzipname = '180228_068';

% averaged bare DNA trace to use for aligning
load('../data/hopping/Hopping/Bare_DNA/averageunzip_processed.mat','avgtrace','distalign','xshift','params')
xshift0 = xshift;
params.dodisplay=0; % avoid displaying graphs for intermediate steps 

opt = struct();
opt.dodisplay=2;

%% load in unzipping trace
load([dirname unzipname])

binw = 2; % smooth over 2 nm
dx = 1; % bins at 1nm spacing
% interpolate to equispaced distance points
nint = 500;
distinterp = linspace(440,1080,nint);
avgtrace0 = interp1(distalign,avgtrace,distinterp);

[binx, fsmooth] = smoothLocalAvg(Ext_raw{1},Fmean_raw{1},binw,dx);
fvals = interp1(binx,fsmooth,distinterp);

% align to the previous averaged curve, to figure out the appropriate shift

ndistalign = 200;
rangealign = [630 870; 630 870];
[distalign,avgtracecur,allfalign,tdel] = getAvgAlignedTrace(distinterp,{avgtrace0,fvals},ndistalign,struct('rangealign',rangealign));

xshift=xshift0-tdel(2)*(distinterp(2)-distinterp(1));

if (opt.dodisplay>1)
    figure
    plot(distinterp,avgtrace0)
    hold all
   % plot(distalign,allfalign{2})
    plot(distinterp+tdel(2)*(distinterp(2)-distinterp(1)),fvals)
    hold off
    xlabel('extension (nm)')
    ylabel('force (pN)')
    legend('avg bare DNA','current')
end

%% Calculate fitted binding energy
params.fix=1; % use parameters determined for avg bare DNA curve
distvals = distinterp-xshift;
[seqshift,fitbindE,params2] = getBindEfitSSParam(distvals,fvals,params);
bindE = fitbindE/params.kT;

seqshift

%% energy landscape
[Evals, Estretch,forces,zplot,nlist] = getElandscapeFromBindE(bindE,params);

% -------------
%% Now process the hopping data
% --------------

%% load in a hopping data set
load([dirname hopname])

Ext = Ext - xshift;
options = struct('dodisplay',opt.dodisplay);
[allpts,meanslope,clustershifts] = clusterParallelLines(Ext,Fmean,options);

%% For each trap position, calculate z(n) and stretching energy for different bp unzipped

[ztrappos,stretchenergytrappos] = fixTrapZEnergies(clustershifts,meanslope,nlist,zplot,forces,Estretch);

%% extract binding energy from equilibrium hopping data
minct=20;
bindEeq = getBindEfromHoppingData(allpts,clustershifts,meanslope,ztrappos,stretchenergytrappos,minct,nlist,params);
if (opt.dodisplay>1)
    figure
    pcolor(bindEeq); shading flat
    caxis([0,4])
    xlabel('bp unzipped')
    ylabel('trap position index')
    title('energy for each bp extracted from eq distrib')
end

%% combine to get average unzipping energies from different trap positions
avgbindE = NaN*zeros(1,size(bindEeq,2));
%avgbindE = bindE;
for nc = 1:size(bindEeq,2)
    ind = find(bindEeq(:,nc)>0 & ~isnan(bindEeq(:,nc)));
    if (~isempty(ind))
        avgbindE(nc) = mean(bindEeq(ind,nc));
    end
end

goodind = find(~isnan(avgbindE));
if (opt.dodisplay)
    figure
    plot(avgbindE)
    xlim([1,nlist(end)])
    xlabel('number unzipped')
    ylabel('binding energy')
end

%%
save([dirname hopname '_processed.mat'])