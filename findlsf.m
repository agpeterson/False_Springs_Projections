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

%% Check accumulated chilling hours below 7.2C base temperature; chilling hours
% should exceed 200/1350 (citation?).
tmin = tmin(1:181); 				% Cut days at June 31.
tmin(tmin > 7.2) = 0;				% Zero all days above 7.2.
chilling_days = -(nansum(tmin));	% Sum with inverse sign.

%% Split on chilling_days, finding latest day in first half of the year.
if chilling_days > 200
	
	% Find latest day of year at or below threshold.
	freeze = max(find(tmin <= threshold));

	% Quality check for NaNs.
	if ~isempty(freeze) & sum(isnan(tmin)) < 20
		lsf = freeze;
	else 
		lsf = NaN;
	end

else

	lsf = NaN;

end 	% If statement.

end 	% Function.