%%=============================================================================
% NAME:   Model_Comparisons.m
% AUTHOR: Alexander Peterson
% DATE:   1 Dec. 2014
%
% DESC:   Calculate 2 standard deviations between historical
% model runs, then find the difference between the multi-model mean
% historical and future. If the difference is greater than 2 standard
% deviations, it's significant. Finally, find pixels where more than 90% of
% models show direction agreement.
% REF:	  None.
% NOTE:	  
%
% IN:     
% OUT:    
% CALL:   
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

% Preallocate mdl_robustness
mdl_robustness = NaN(N_LAT,N_LON,'single');

% Iterate over lat/lon and check model robustness. if absolute value of multi-model mean change
% is greater or equal to 2*stnd;
for x=1:N_LON
	for y=1:N_LAT

		tmp = squeeze(fut_diff(y,x,:));
		sign_sum = nansum(sign(tmp));
		sign_pcnt = sign_sum / N_MDL;

		if abs(fut_diff_mean(y,x)) >= hst_two_stnd(y,x) && ...
				abs(sign_pcnt) >= 0.9
			mdl_robustness(y,x) = 1;
		elseif abs(fut_diff_mean(y,x)) < hst_two_stnd(y,x) || ...
				abs(sign_pcnt) < 0.9
			mdl_robustness(y,x) = 0;
		else
			mdl_robustness(y,x) = NaN;
		end
	end
end


end