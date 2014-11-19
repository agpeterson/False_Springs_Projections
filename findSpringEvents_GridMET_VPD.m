%%=============================================================================
% NAME:   findSpringEvents.m
% AUTHOR: John Abatzoglou, Alexander Peterson, Katherine Hegewisch
% DATE:   1 Nov. 2014
%
% DESC:   The function iterates in parallel over the GridMET dataset
%		  and returns VPD-derived green-up dates as subsets of CONUS.
% REF:	  None.
% NOTE:	  Inputs can be class-type single or double.
%
% IN:     tmin - Daily minimum temperature in units K in row/column vector.
%		  tmax - Daily maximim temperature in units K in row/column vector.
%		  sph - Daily specific humidity.
%		  pres - Elevation-derived atmospheric pressure.
%		  photoperiod - Photoperiod for all calendar days in row/column vector.
%		  N_DAYS - Number of days as a scalar value.
%		  N_YRS - Number of years as a scalar value.
%		  N_LAT - Number of points of latitude to iterate over as scalar value.
%		  SCN - Scenario number between 1 and 3 indicating outer script 
%				progress; used to calculate or use gsi_max.
% OUT:    gu_sub - Subset of green-up dates.
% CALL:   findLSF.m; calcGSI.m; findGreenup.m
%==============================================================================


function [gu_sub] = findSpringEventsVPD(tmin,tmax,sph,pres,photoperiod,...
										N_DAYS,N_YRS,N_LAT,SCN);


%%=============================================================================
% Constants, initialized variables, and data control.
%==============================================================================
% Preallocate arrays for LSF and GSI geographical subsets.
gsi_sub = NaN(N_YRS,N_LAT,'single');
gsi_max = NaN(N_LAT,'single');

% Convert rom Kelvin to C and cast inputs as doubles.
tmin = double(tmin - 273.15);
tmax = double(tmax - 273.15);
sph = double(sph);


%%=============================================================================
% Main loop.
%==============================================================================
% Iterate over lat.
for lat=1:N_LAT

	% Find non-empty latitudes.
	f = find(~isnan(tmin(:,:,lat)));
	
	% If tmin is not empty, iterate over years.
	if ~isempty(f)
		
		for yr=1:N_YRS

			% Call dew point function.
            tdew = calcDewPoint(sph(:,yr,lat),pres(lat));

            % Call VPD function.
            vpd = calcVPD(tmax(:,yr,lat),tmin(:,yr,lat),tdew,pres(lat));

			% Call calcGSI.m function.
			gsi_raw(:,yr) = calcGSI(tmin(:,yr,lat),photoperiod(:,lat),vpd);
	    	
	    end 	% yr; N_YRS.
		
		% Call findGreenup.m function and concatenate to sub-region.
		[gu_sub(:,lat),gsi_max(lat)]= findGreenup(gsi_raw,gsi_max(lat),...
												   N_DAYS,N_YRS,SCN);
	
	end 	% If statement.
end         % lat; lat_subset.





end 		% Function.