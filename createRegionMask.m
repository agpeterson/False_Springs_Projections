%%=============================================================================
% NAME:   findLatLonIndices.m
% AUTHOR: Alexander Peterson
% DATE:   4 Jan. 2015
% DESC:   The function finds regional indices from given lat and lon grids,
%			which can be created using the meshgrid() function.
% IN:     rgn_coords - an X by 2 array of coordinate pairs.
%		  lon_grid -
%		  lat_grid -
% OUT:    
% CALL:   
%==============================================================================

function [rgn_mask] = createRegionMask(rgn_coords,lon_grid,lat_grid)

%%=============================================================================
% Iterate over regional coordinates, finding matches within error tolerance.
disp('Iterating over regional coordinates and creating index...')
for i=1:length(rgn_coords)
	tmp = find(rgn_coords(i,1) >= lon_grid(:,:) - 0.0001 & ...
			   rgn_coords(i,1) <= lon_grid(:,:) + 0.0001 & ...
			   rgn_coords(i,2) >= lat_grid(:,:) - 0.0001 & ...
			   rgn_coords(i,2) <= lat_grid(:,:) + 0.0001);
	
	% Continue if find returns empty, otherwise store index value.
	if isempty(tmp)
		continue
	else
		rgn_index(i,:) = tmp;
	end % if statement

end 	% i; 1:length(rgn_coords)


%%=============================================================================
disp('Creating lat x lon regional grid...')

% Create variables to be used in creating regional grid of 0s and 1s.
n_lat = size(lon_grid,1);
n_lon = size(lon_grid,2);

% Create gridded cell array to be indexed into.
rgn_mask = num2cell(zeros(n_lat,n_lon));

% Iterate over length of regional index and replace corresponding 0s with 1s
% to create regional mask to be returned from function.
for i=1:length(rgn_index)
	rgn_mask{rgn_index(i)} = 1;
end

% Convert from cell to numeric array.
rgn_mask = cell2mat(rgn_mask);
disp('Function complete.')


end 	% function