function [binx,fsmooth] = smoothLocalAvg(distvals,fvals,binw,dx,binrange)
%% Smooth force-extension data over distances 
% average over bins of width binw, spaced at intervals dx
if (nargin>4)
    binx = binrange(1):dx:binrange(2);
else
    binmin = min(distvals)+binw+dx;
    binmax = max(distvals)-binw-dx;
    binx = binmin:dx:binmax;
end
nb = length(binx);
for bc = 1:nb
    %%
    m1 = binx(bc)-binw/2;
    m2 = m1+binw;
    ind = find(distvals>m1 & distvals<=m2);
    fsmooth(bc) = mean(fvals(ind));
end

% fill in NaNs by interpolation
goodind = find(~isnan(fsmooth));
badind = find(isnan(fsmooth));

fsmooth(badind) = interp1(binx(goodind),fsmooth(goodind),binx(badind));