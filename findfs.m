%%=============================================================================
% NAME:   findfs.m
% AUTHOR: Alexander Peterson
% DATE:   21 Sept. 2014
% DESCR:  This function finds false springs.
% IN:     GSI, LSF, and a day lag.
% OUT:    FS
% CALLS:  N/A
%==============================================================================

function [fs,fsei] = findfs(gsi,lsf,LAG)

% Create variable for number of years.
n_yrs = size(gsi,1);

% Iterate over years to classify false springs as binary, where 1 is LSF
% occurring LAG or more days post-greenup.
for i=1:n_yrs

	if isnan(lsf(i)) || isnan(gsi(i))
		fs(i) = NaN;
	elseif lsf(i) >= LAG + gsi(i)
		fs(i) = 1;
	else
		fs(i) = 0;
	end 	% If statement.

end 		%i; n_yrs.

% Sum number of false springs.
fs_sum = nansum(fs(:));

% Find exposure index.
fsei = (fs_sum / sum(~isnan(lsf)))*100;

end 		% Function.
