%%=============================================================================
% NAME:   findLSF.m
% AUTHOR: John Abatzoglou, Alexander Peterson, Katherine Hegewisch
% DATE:   1 Nov. 2014
%
% DESC:   This function finds the last spring freeze date for one year of daily
%		  minimum temperature.
% REF:	  Richardson et al. (1975)
% NOTE:	  Inputs can be class-type single or double.
%
% IN:     tmin - Daily minimum temperature in units C in row/column vector.
%		  THRESHOLD - Threshold constant defining spring freeze as scalar value.
% OUT:    lsf - Last calendar date meeting or exceeding threshold.
% CALL:   N/A
%==============================================================================


function [lsf] = findlsf(tmin,THRESHOLD)


%%=============================================================================
% Find latest day of year at or below threshold.
%==============================================================================
lsf = max(find(tmin <= THRESHOLD));

% If freeze is not empty, set lsf equal to freeze value. If empty, set to 0.
if isempty(lsf)
	lsf = single(0);
end



end 	% Function.