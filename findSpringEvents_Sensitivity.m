%%=============================================================================
%
% NAME:   findSpringEvents.m
% AUTHOR: Alexander Peterson
% DATE:   1 Nov. 2014
% DESC:   The function iterates in parallel over the MACAv2-METDATA dataset
%		  and returns last spring freeze and green-up dates as subsets of
%		  CONUS.
% IN:     tmin - Daily minimum temperature in units K in row/column vector.
%		  gsi_max - Max GSI values pre-calculated or empty array.
%		  photoperiod - Photoperiod for all calendar days in row/column vector.
%		  vpd - Vapor pressure deficit in row/column vector.
%		  N_DAYS - Number of days as a scalar value.
%		  N_YRS - Number of years as a scalar value.
%		  N_LAT - Number of points of latitude to iterate over as scalar value.
%		  HST_INDEX - Indices of historical years.
%		  SCN - Scenario number between 1 and 3 indicating outer script 
%				progress; used to calculate or use gsi_max.
% OUT:    lsf_sub - Subset of last spring freeze dates.
%		  gu_sub - Subset of green-up dates.
%		  gsi_max - Max GSI values to be used in normalization.
% CALL:   findLSF.m; calcGSI.m; findGreenup.m
%
%==============================================================================


function [lsf_sub,gu_sub] = findSpringEvents(tmin,delta,...
	photoperiod,vpd,N_DAYS,N_YRS,N_LAT,SCN);

% Uncomment if using VPD in GSI.
% function [gu_sub] = findSpringEventsVPD(tmin,tmax,sph,pres,photoperiod,...
%	N_DAYS,N_YRS,N_LAT,SCN);


%%=============================================================================
% Constants, initialized variables, and data control.
%==============================================================================
% Preallocate arrays for LSF and GSI geographical subsets.
lsf_sub = NaN(N_YRS,N_LAT,'single');
gsi_raw = NaN(N_DAYS,N_YRS,'single');
gu_sub = NaN(N_YRS,N_LAT,'single');
gsi_max = NaN(N_LAT,1,'single');


% Convert from Kelvin to C and cast inputs as doubles.
tmin = double(tmin - 273.15);
delta = double(delta);

% Uncomment if using VPD in GSI.
% tmax = double(tmax - 273.15);
% sph = double(sph);

% Set LSF threshold.
THRESHOLD = single(-2.2);

% Modify temperatures for sensitivity experiment.
for lat=1:N_LAT
	tmin_mod(:,:,lat) = tmin(:,:,lat) + delta(lat,:);
end


%%=============================================================================
% Main loop.
%==============================================================================
% Iterate over lat.
for lat=1:N_LAT

	% Find non-empty latitudes.
	f = find(~isnan(tmin_mod(:,:,lat)));
	
	% If tmin is not empty, iterate over years.
	if ~isempty(f)
		
		for yr=1:N_YRS

			% Call LSF function.
			lsf_sub(yr,lat) = findLSF(tmin_mod(:,yr,lat),THRESHOLD);

			% Call dew point function - uncomment if using VPD in GSI.
            % tdew = calcDewPoint(sph(:,yr,lat),pres(lat));

            % Call VPD function - uncomment if using VPD in GSI.
            % vpd = calcVPD(tmax(:,yr,lat),tmin(:,yr,lat),tdew,pres(lat));

			% Call calcGSI.m function.
			gsi_raw(:,yr) = calcGSI(tmin_mod(:,yr,lat),photoperiod(:,lat),vpd);
	    	
	    end 	% yr; N_YRS.
		
		% Call findGreenup.m function and concatenate to sub-region.
		[gu_sub(:,lat),gsi_max(lat)]= findGreenup(gsi_raw,gsi_max(lat),...
												   N_DAYS,N_YRS,SCN);
	
	end 	% If statement.
end         % lat; lat_subset.


end 		% Function.