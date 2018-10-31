function [distalign,avgtrace,allfalign,tdel] = getAvgAlignedTrace(distinterp,allfinterp,ndistalign,options)
% align all curves to the first one
% assumes all curves have been interpolated to the same distinterp values
% returns distances and aligned force curves
% also returns the average trace

% default ooptions
opt = struct();
% range of extension values to align, for each trace
opt.rangealign = zeros(length(allfinterp),2);
opt.rangealign(:,1) = 630; 
opt.rangealign(:,2) = 980;

%%
if (nargin>3)
    opt = copyStruct(options,opt);
end
%%

fc1 = 1;

%plot(distinterp,allfinterp{fc1},'k','LineWidth',2)    
%hold all

for fc2 = 1:length(allfinterp)
    
    ind = find(distinterp>opt.rangealign(fc2,1) & distinterp<opt.rangealign(fc2,2));
    % ind = find(distinterp>630 & distinterp<800);
    s1 = allfinterp{fc1}(ind);
    s2 = allfinterp{fc2}(ind);
    tvals = distinterp(ind);
    
    %
    s1 = s1 - mean(s1);
    s2 = s2 - mean(s2);
    
    % maximum shift allowed
    maxlag = 30;
    
    tdel(fc2) = finddelay(s2,s1,maxlag);        
end

%% create aligned versions of all the traces
mint = min(tdel); maxt = max(tdel);

% aligned traces will be evaluated at these distances
distalign = linspace(distinterp(maxt+1),distinterp(end-abs(mint)),ndistalign);
scl = 1;

avgtrace = zeros(size(distalign));
for fc2 = 1:length(allfinterp)
     if (tdel(fc2)>0)      
        allfalign{fc2} = interp1(distinterp(tdel(fc2)+1:end),allfinterp{fc2}(1:end-tdel(fc2))*scl,distalign);
     else       
        allfalign{fc2} = interp1(distinterp(1:end-abs(tdel(fc2))),allfinterp{fc2}(abs(tdel(fc2))+1:end)*scl,distalign);
     end
    
     avgtrace = avgtrace + allfalign{fc2};
end
avgtrace = avgtrace/length(allfalign);

end