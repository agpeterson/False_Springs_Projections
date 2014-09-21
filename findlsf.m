%%=============================================================================
% NAME:   findlsfl.m
% AUTHOR: Alexander Peterson
% DATE:	  12 Sept. 2014
% DESCR:  This function finds the last spring freeze date for one year of daily
% 		  minimum temperature data.
% IN:     Daily minimum temperature data with dimensions [DAYS].
% OUT:    Julian day of last freeze with dimensions [LSF].
% CALLS:  No external functions.
%==============================================================================

function [lsf] = findlsf(tmin,threshold)

% Restrict days at June 31, Julian Date 181 to avoid running into fall freezes.
tmin = tmin(1:181);

% Find latest day of year at or below threshold.
freeze = max(find(tmin <= threshold));

% If freeze is not empty, set lsf equal to freeze value. If empty, set to 0.
if ~isempty(freeze)
	lsf = freeze;
else
	lsf = 0;
end

% Quality check for missing values; average over temperature - if NaN, set lsf
% equal to NaN.
tmin_mean = mean(tmin);
if isnan(tmin_mean)
	lsf = NaN;
end

end 	% Function.
