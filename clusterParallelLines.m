function [finalptdata, meanslope,clustershifts, linecoeff] = clusterParallelLines(Ext,Fmean,options)
% cluster data consisting of X= Ext values and Y = Fmean values
% into groups that approximate parallel lines
% returns cell array of x,y data (finalptdata)
% meanslope = average slope of all the parallel lines
% clustershifts = estimate of shift for each cluster
% linecoeff = coefficients of individual linear fits

% set default parameters

% if dodisplay = 0, no plots
% if dodisplay = 1, final plot only
% if dodisplay = 2, make intermediate plots
opt.dodisplay = 1; 
opt.firstind = 1:2000; % indices within first group, used to get slope estimate
opt.minpeakheight = 100; % minimum number of data points defining histogram peak
opt.zcutoff = 5; % number of std dev used to get rid of outliers in clusters

% copy over input parameters
if (nargin>2)
    opt = copyStruct(options,opt);
end

%% Fit early data points to get slope approximation
firstind = opt.firstind;
cfit = polyfit(Ext(firstind),Fmean(firstind),1);
slopeapprox = cfit(1);
% get perpendicular line 
p0 = [mean(Ext(firstind)) mean(Fmean(firstind))]; % center point
slopeperp = -1/slopeapprox;
bperp = p0(2) - slopeperp*p0(1);

if (opt.dodisplay>1)
    figure
    plot(Ext,Fmean,'.')
    hold all
    plot(Ext(firstind),Fmean(firstind),'.')
    extlist = linspace(min(Ext(firstind))-10,max(Ext(firstind))+10);
    plot(extlist,cfit(1)*extlist+cfit(2))
    %plot(extlist,slopeperp*extlist+bperp)
    hold off
end
%% project all data points onto the perpendicular line
pts = [Ext' Fmean'];
ptshift = bsxfun(@plus, pts,-p0);
v = [1;slopeperp]; % vector giving line direction
ptproj = ptshift*v/norm(v)^2;

%% histogram the projected points and find peaks
nbin = 200;
[cts,binx] = hist(ptproj,nbin);
dx = binx(2)-binx(1);
% which points go in which bin?
edges = binx(1)-dx/2:dx:binx(end)+dx/2;
binids = discretize(ptproj,edges);

[pks,locs] = findpeaks(cts);
ind = find(pks>opt.minpeakheight);
pks = pks(ind); locs = locs(ind);

if (opt.dodisplay>1)
    figure
    plot(binx,cts,'.-',binx(locs),pks,'*')
end

%% classify all data points by pulling out those within 
% half of minimal-interpeak-distance of each peak
maxpeakwidth = min(diff(binx(locs)))/2;
%plot(Ext,Fmean,'k.')
if (opt.dodisplay>1)
    plot(1:length(ptproj),ptproj,'k.')
end
ptgroups = {};
for lc = 1:length(locs)
    closebins = find(abs(binx-binx(locs(lc))) < maxpeakwidth);
    ptgroups{lc} = find(ismember(binids,closebins));
    if (opt.dodisplay>1)
    hold all
    plot(ptgroups{lc},ptproj(ptgroups{lc}),'.')
%    plot(Ext(ptgroups{lc}),Fmean(ptgroups{lc}),'.')
    end
end
if (opt.dodisplay>1)
    hold off
end
ptgroups0 = ptgroups;

%% get rid of outlier points
for gc  = 1:length(ptgroups)
    curpts = ptproj(ptgroups0{gc});
    s = std(curpts);
    mn = mean(curpts);
    
    sclpts = (curpts-mn)/s;
    zscorepts = zscore(sclpts);
    ptgroups{gc} = ptgroups0{gc}(abs(zscorepts)<opt.zcutoff);
end

%% save grouped point data
for gc = 1:length(ptgroups)
        finalptdata{gc} = pts(ptgroups{gc},:);       
end


%% fit a line to each cluster of points
linecoeff = zeros(length(finalptdata),2);
for gc = 1:length(finalptdata)
    pts = finalptdata{gc};
    linecoeff(gc,:) = polyfit(pts(:,1),pts(:,2),1);
end

meanslope = -mean(linecoeff(:,1));

if (opt.dodisplay>1)
    figure
    plot(Ext,Fmean,'k.')
    hold all
    for gc = 1:length(ptgroups)     
        plot(finalptdata{gc}(:,1),finalptdata{gc}(:,2),'.')
        extlist = linspace(min(finalptdata{gc}(:,1)) - 5,max(finalptdata{gc}(:,1)) + 5)
        plot(extlist,linecoeff(gc,1)*extlist+linecoeff(gc,2),'k');
        hold all
    end
    hold off
end

%% fit a shift for each cluster, assuming constant slope
clustershifts = zeros(length(finalptdata),1);
for gc = 1:length(finalptdata)
    pts = finalptdata{gc};
    clustershifts(gc)= mean((pts(:,2) + meanslope*pts(:,1))/meanslope);
end

if (opt.dodisplay>0)
    figure
    cmat = jet(length(finalptdata))
    plot(Ext,Fmean,'k.')
    hold all
    for gc = 1:length(finalptdata)     
        plot(finalptdata{gc}(:,1),finalptdata{gc}(:,2),'.','Color',cmat(gc,:))
        extlist = linspace(min(finalptdata{gc}(:,1)) - 5,max(finalptdata{gc}(:,1)) + 5)
        plot(extlist,meanslope*(clustershifts(gc) - extlist),'k');
        hold all
    end
    hold off
end

