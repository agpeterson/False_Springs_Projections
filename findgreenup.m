%%=============================================================================
% NAME:   findgreenup.m
% AUTHOR: Alexander Peterson
% DATE:   16 Sept. 2014
% DESCR:  This function processed raw daily GSI data.
% IN:     Raw GSI values with array dimensions [DAYS,YEARS].
% OUT:    
% CALLS:  movingaverage.m
%==============================================================================

function [gsi] = findgreenup(gsi_raw)

% Create variables for number of years and latitude, and smoothing window.
n_day = size(gsi_raw,1);
n_yrs = size(gsi_raw,2);
window = 10;

% Reshape raw GSI data to continuous time-series.
gsi_continuous = reshape(gsi_raw,n_day*n_yrs,1);

% Call moving average function across continuous time-series.
gsi_mvg_avg = movingaverage(gsi_continuous,window);

% Find 95th percentile of the moving average for normalization.
gsi_max = prctile(gsi_mvg_avg,95,2);

% Reshape continuous time-series back to days, years.
gsi_avg = reshape(gsi_mvg_avg,n_day,n_yrs);

% Divide average GSI values by 95th percentile maximum to normalize.
gsi_norm = gsi_avg./gsi_max;

% Iterate over years to find green-up.
for i=1:n_yrs

	% Find earliest day of year where GSI meets green-up threshold.
	greenup = min(find(gsi_norm(:,i) >= 0.5));
    
    if ~isempty(greenup)
        gsi(:,i) = greenup;
    else
        gsi(:,i) = NaN;
    end 	% If statement.

end 		% i; n_yrs.

end 		% Function loop.