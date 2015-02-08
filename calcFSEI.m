%%===========================================================================
%
% NAME:   calcFSEI.m
% AUTHOR: Alexander Peterson
% DATE:   22 Jan. 2015
% DESC:   Classifies false springs and derives FSEI.
%
%============================================================================

function [fs,fsei] = calcFSEI(gu,lsf,DAY_LAG)


%%============================================================================
% Initial variables and constants.
%-----------------------------------------------------------------------------

N_YRS = size(gu,1);
N_LAT = size(gu,2);
N_LON = size(gu,3);


%%============================================================================
% Classify false springs.
%-----------------------------------------------------------------------------

% Preallocate FS variable.
disp('Preallocating FS variable...')
fs = NaN(N_YRS,N_LAT,N_LON);

% Iterate over lon, lat, and years.
disp('Classifying false springs...')
for lon=1:N_LON
	for lat=1:N_LAT
		for yrs=1:N_YRS
			if isnan(lsf(yrs,lat,lon)) || isnan(gu(yrs,lat,lon))
				fs(yrs,lat,lon) = NaN;
			elseif lsf(yrs,lat,lon) >= (DAY_LAG + gu(yrs,lat,lon))
				fs(yrs,lat,lon) = 1;
			else
				fs(yrs,lat,lon) = 0;
			end 	% if statement
		end 	% yrs; 1:N_YRS
	end 	% lat; 1:N_LAT
end 	% lon; 1:N_LON


%-----------------------------------------------------------------------------
% Classify false springs.
%-----------------------------------------------------------------------------

% Preallocate FSEI variables.
disp('Preallocating FSEI variables...')
fs_sum = NaN(N_LAT,N_LON);
gu_sum = NaN(N_LAT,N_LON);
fsei = NaN(N_LAT,N_LON);

% Derive FSEI.
disp('Deriving FSEI...')
fs_sum(:,:) = nansum(fs(:,:,:),1);	% Sums ones.
gu_sum(:,:) = sum(~isnan(gu(:,:,:)),1); % Sums # of elements.
fsei(:,:) = (fs_sum(:,:) ./ gu_sum(:,:)) * 100;


end 	% Function.