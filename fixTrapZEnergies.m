function [ztrappos,stretchenergytrappos] = fixTrapZEnergies(clustershifts,meanslope,nlist,zplot,forces,Estretch)
% given fixed trap positions (in clustershifts)
% calculate the z extension associated with each length unzipped
% and also the stretching energy associated with each length unzipped
% NOTE: clustershifts actually contain trap position - bead radii
% meanslope = trap stiffness (slope of clusters)
% nlist, zplot = lists of unzipped basepairs and z extension over which
% the landscapes are calculated
% forces = landscape of forces associated with each n, z combination
% Estretch = energy landscape for stretching energy only

[N,Z] = meshgrid(nlist,zplot);
for sc = 1:length(clustershifts);
    shift = clustershifts(sc);
    
    for nc = 1:length(nlist) % for each number of bp unzipped
        % find where chain force = trapforce
        forcediff = forces(nc,:) - meanslope*(shift-zplot);
        ind = find(~isnan(forcediff));
        ztrappos(sc,nc) = interp1(forcediff(ind),zplot(ind),0,'pchip');        
    end
    % interpolate energies along this curve with fixed trap position    
    stretchenergytrappos(sc,:) = interp2(N,Z,Estretch',nlist,ztrappos(sc,:));
end

