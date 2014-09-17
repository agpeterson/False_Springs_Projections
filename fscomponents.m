%%=============================================================================
% NAME:   fscomponents.m
% AUTHOR: Alexander Peterson
% DATE:   12 Sept. 2014
% DESCR:  This function iterates in parallel over the MACA dataset to determine
%		  last spring freeze and GSI.
% IN:     Arrays t_var, day_subset, vpd_subset, n_yrs, lat_subset.
% OUT:    Arrays lsf_sub and gsi_sub.
% CALLS:  findlsf.m and calcgsi.m
%==============================================================================

function [lsf_sub,gsi_sub] = fscomponents(tmin,day_subset,vpd_subset);

% Create variables for number of years and latitude.
n_day = size(tmin,1);
n_yrs = size(tmin,2);
n_lat = size(tmin,3);

% Preallocate variables for subregional LSF and GSI.
lsf_sub = NaN(n_yrs,n_lat);
gsi_sub = NaN(n_yrs,n_lat);
gsi_raw = NaN(n_day,n_yrs,n_lat);

% Iterate over latitude and years, calling other functions.
for i=1:n_lat

	% Create temporary variable for temperature data and find all NaNs.
	t_var = tmin(:,:,i);
	f = find(~isnan(t_var));
	
	% If t_var is not empty, iterate over years to process LSF and GSI. 
	if(~isempty(f))
		
		for j=1:n_yrs
		
			% Call findlsf and calcgsi function.
			lsf_sub(j,i) = findlsf(t_var(:,j),-2.2);	% 270.95 Kelvin.
			gsi_raw(:,j,i) = calcgsi(t_var(:,j),...
									 day_subset(:,i),...
									 vpd_subset(:,i));
	    
	    end 	% j; n_yrs.
		
		% Call findgsi function and concatenate to sub-region.
		t_gsi = gsi_raw(:,:,i);		% Create temporary variable for GSI.
		gsi_sub(:,i) = findgsi(t_gsi);

	end 	% If statement.

end         % i; lat_subset.

end 		% Function loop.
