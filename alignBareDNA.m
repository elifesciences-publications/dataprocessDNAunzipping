%% Use this code to align several unzipping traces by shifting along the extension axis

%% Data files for bare DNA
dirname = '../data/hopping/Hopping/Bare_DNA/'; % directory where data files are stored
% file names
fnames = {'180305_009.mat','180305_022.mat','180305_032.mat','180305_054.mat','180305_062.mat','180306_078.mat','180312_050.mat'};

%% Load in all unzipping curves
data=cell(1,length(fnames));

binw = 2; % smooth over 2 nm
dx = 1; % bins at 1nm spacing

% interpolate to equispaced distance points
nint = 500;
distinterp = linspace(440,1080,nint);
alltvals = {}; alldistvals = {}; allfvals = {};
for fc = 1:length(fnames)        
    
    load([dirname fnames{fc}]);
    alldistvals{fc} = Ext_raw{1};
    allfvals{fc} = Fmean_raw{1};
    
    % look at initial extension only
    [a,b] = max(alldistvals{fc});
    %ind = b:length(alldistvals{fc});        
    ind=1:b;
    %ind = 1:length(alldistvals{fc});
    
    [allbinx{fc} allfsmooth{fc}] = smoothLocalAvg(alldistvals{fc}(ind),allfvals{fc}(ind),binw,dx);
    allfinterp{fc} = interp1(allbinx{fc},allfsmooth{fc},distinterp);
    distvals = alldistvals{fc}; fvals = allfvals{fc};
    
    plot(distvals(ind),fvals(ind),'.')
    hold all
    drawnow
end
hold off

%% Align all the curves
ndistalign = 200;
[distalign,avgtrace,allfalign,tdel] = getAvgAlignedTrace(distinterp,allfinterp,ndistalign);
for fc = 1:length(allfalign)
    plot(distalign,allfalign{fc})
    hold all
end
hold off

%% work with the average trace
fvals = avgtrace; 
distinterp = distalign;

%% find appropriate shift along x axis, so that the repeats in
% the calculated binding energy are separated by 197 bp
% fit ss DNA parameters (bss, ssdel) as you go

guessshift=-10;

params = struct();
params.lp=35.4; 
params.kT=4.1; 
params.K=1020; 
params.Lds = 924*0.34;
params.Nss = 434; % length of hairpin
% total number of ss bp when unzipped
params.Nsstot = (2*params.Nss+4);
params.Kss = 1000;
% expected start and end points of the sequence repeats
params.start1 = 30; 
params.end1 = 78;
params.start2 = 30+197;
params.end2 = params.start2 + params.end1-params.start1;

options = struct('display','iter');
wantseqshift = 197; % desired shift in bp for the binding profile

% initial guesses for persistence length bss and 
% contour length per bp ssdel
params.bss = 1;
params.ssdel=0.6;
params.Lss = params.ssdel*params.Nsstot;
params.fix = 0;
params.dodisplay=0;

% xshift is the number of nm to shift the extension to the left
xshift = bisectionsolve(@(s) getBindEfitSSParam(distinterp-s,fvals,params)-wantseqshift,[guessshift-5,guessshift+5],options)

%% get the fitted parameters for this shift
distvals = distinterp-xshift; 

params = struct();
params.lp=35.4; 
params.kT=4.1; 
params.K=1020; 
params.Lds = 924*0.34;
params.Nss = 434; % length of hairpin
% total number of ss bp when unzipped
params.Nsstot = (2*params.Nss+4);
params.Kss = 1000;

% expected start and end points of the sequence repeats
params.start1 = 20; 
params.end1 = 78;
params.start2 = 20+197;
params.end2 = params.start2 + params.end1-params.start1;
params.dodisplay=2;

[seqshift,fitbindE,params] = getBindEfitSSParam(distvals,fvals,params);
seqshift
%% Save aligned trace, xshift, and fitted ssDNA parameters
save('../test/averageunzip_processed.mat','avgtrace','distalign','xshift','params')
