%%=============================================================================
% NAME:   parallel_function.m
% AUTHOR: Alexander Peterson
% DATE:   12 Sept. 2014
% DESCR:  This function iterates in parallel over the MACA dataset to determine
%		  last spring freeze and GSI.
% IN:     Arrays tmin, day_subset, vpd_subset, n_yrs, lat_subset.
% OUT:    Arrays lsf_sub and gsi_sub.
% CALLS:  findlsf.m and calcgsi.m
%==============================================================================

function [lsf_sub,gsi_sub] = parallel_function(tmin,day_subset,vpd_subset,...
											   n_yrs,lat_subset)

% Iterating over lat subset and years, calling other functions.
for l=1:length(lat_subset)
	for y=1:n_yrs
		% Create temporary variable to hold lat and year subset.
		t_var = tmin(:,y,l);

		% Call findlsf and calcgsi function.
		lsf_sub(y,l) = findlsf(t_var,-2.2);
		gsi_sub(y,l) = calcgsi(t_var,day_subset(:,l),vpd_subset(:,l));
	end     % y; n_yrs.
end         % l; lat_subset.

end 		% Function loop.