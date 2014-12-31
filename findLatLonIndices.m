%%=============================================================================
% NAME:   findLatLonIndices.m
% AUTHOR: Alexander Peterson
% DATE:   30 Dec. 2014
%
% DESC:   The function finds subset indices from given lat and lon.
% REF:	  None.
% NOTE:	  Inputs can be class-type single or double.
%
% IN:     region - Regional lat and lons to be 
% OUT:    
% CALL:   
%==============================================================================

function [lat_index,lon_index] = findLatLonIndices(region,lat,lon)

% Find unique lat/lon points from grid.
lon_new = unique(region(:,1));
lat_new = unique(region(:,2));

% Iterate over all lats, finding lats corresponding to regional lats. Create
% logical index where lats match.
for i=1:length(lat)
	lat_compare = find(lat_new >= (lat(i)-.01) & lat_new <= (lat(i)+.01));
	if isempty(lat_compare)
		lat_index(i) = 0;
	else
		lat_index(i) = 1;
	end
end

% Same as above, but for lons.
for i=1:length(lon)
	lon_compare = find(lon_new >= (lon(i)-.01) & lon_new <= (lon(i)+.01));
	if isempty(lon_compare)
		lon_index(i) = 0;
	else
		lon_index(i) = 1;
	end
end

% Create lat and lon indices to be used in subsetting.
lat_index = find(lat_index == 1);
lon_index = find(lon_index == 1);

end