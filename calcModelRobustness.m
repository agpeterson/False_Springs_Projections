%%=============================================================================
% NAME:   caclModelRobustness.m
% AUTHOR: Alexander Peterson
% DATE:   17 Jan. 2015
% DESC:   Determines model significance across a gridded geographic dataset,
%			return 1 where the multi-model mean is greater than 2 standard
%			devations of historical model values and where more than 90% of
%			models agree on direction of change (negative/positive).
% REF:	  IPCC 2013.
% IN:     hst_mean - [lat,lon,models]
%		  fut_diff - [lat,lon,models]
% OUT:    mdl_robustness - [lat,lon]
% CALL:   N/A
%==============================================================================


function [mdl_robustness] = calcModelRobustness(hst_mean,fut_diff)


%%=============================================================================
% Set variables and constants.
%==============================================================================

N_LAT = size(hst_mean,1);
N_LON = size(hst_mean,2);
N_MDL = size(hst_mean,3);


%%=============================================================================
% Calculate standard deviation among models for historical 1950-2005 period.
%==============================================================================

% Preallocate.
hst_stnd = NaN(N_LAT,N_LON,'single');

% Iterate over lat/lon and call nanstd() function.
for x=1:N_LON
	for y=1:N_LAT
		hst_stnd(y,x) = nanstd(squeeze(hst_mean(y,x,:)));
	end 	% y; 1:N_LAT
end 	% x; 1:N_LON

% Create 2*stnd variable.
hst_two_stnd = hst_stnd * 2;


%%=============================================================================
% Create multi-model historical and futures difference.
% NOTE: Per-model differences previously calculated.
%==============================================================================

% Create multi-model mean.
fut_diff_mean = squeeze(nanmean(fut_diff,3));


%%=============================================================================
% Determine model robustness per pixel by finding where the multi-model mean
% future change is greater than two standard deviations and where more than
% 90% of models agree on direction of change.
%==============================================================================

% Preallocate.
mdl_robustness = NaN(N_LAT,N_LON,'single');

% Iterate over lat/lon; create temporary variable to hold all model change
% values and sum, then calculate the percent. Check if absolute value of the
% multi-model mean change is greater than historical two standard deviations
% and if more than 90% of models agree on direction change.
for x=1:N_LON
	for y=1:N_LAT

		tmp = squeeze(fut_diff(y,x,:));
		sign_sum = nansum(sign(tmp));
		sign_pcnt = sign_sum / N_MDL;

		if abs(fut_diff_mean(y,x)) >= hst_two_stnd(y,x) && abs(sign_pcnt) >= 0.9
			mdl_robustness(y,x) = 1;
		elseif abs(fut_diff_mean(y,x)) < hst_two_stnd(y,x) || abs(sign_pcnt) < 0.9
			mdl_robustness(y,x) = 0;
		else
			mdl_robustness(y,x) = NaN;
		end 	% if statement
	
	end 	% y; 1:N_LAT
end 	% x; 1:N_LON


end 	% Function