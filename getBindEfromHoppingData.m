function bindEeq = getBindEfromHoppingData(allpts,clustershifts,meanslope,ztrappos,stretchenergytrappos,minct,nlist,params)
% given a force-extension data for many observations with a fixed trap
% position,
% calculate the binding energy associated with each bp that shows up with
% sufficiently high probability
% repeat for all trap positions
% allpts = cell array, giving extensions and forces for each trap position
% clustershifts = trap positions minus bead radii 
% meanslope = trap stiffness
% ztrappos = for each trap position, z extension associated with each bp
% unzipped
% stretchenergytrappos = for each trap position stretching energy
% associated with each bp unzipped
% minct = minimum number of observations required to calculate energy for
% this bp for a given trap position
% nlist = full list of number unzipped for which energies can be calculated
% params = mechanical parameters for ss and ds DNA


bindEeq = zeros(length(clustershifts),length(nlist));

for sc = 1:size(bindEeq,1)
    % get the observed probability distribution
    zvals = allpts{sc}(:,1)';
    fvals = allpts{sc}(:,2)';
   
    [freq,binmid,cts] = histogramNunzip(zvals,fvals,nlist,params);
    
    
    % observed energies from prob distrib
    Eobs = -log(freq);
    ind = find(cts>minct);
    ind = ind(1):ind(end);
    Eobs = Eobs(ind);
    nfound = nlist(ind);
    
    % subtract off trap energy and chain stretching energy
    Ebindonly = Eobs - stretchenergytrappos(sc,ind) -  0.5*meanslope*(clustershifts(sc)-ztrappos(sc,ind)).^2/params.kT;
   
    bindEeq(sc,ind(2:end)) = diff(Ebindonly);
end
