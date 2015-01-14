%%=============================================================================
% NAME:   findGreenup.m
% AUTHOR: John Abatzoglou, Alexander Peterson, Katherine Hegewisch
% DATE:   1 Nov. 2014
%
% DESC:   This function processes raw daily GSI data, as calculated by the
%		  calcGSI.m function, and outputs the first calendar day where the 21
%	 	  day moving average meets or exceeds 0.5, classified as 'green-up'.
%		  This function was written to iterate over the MACAv2-METDATA dataset,
%		  and uses the 1950-2005 historical to normalize GSI at all points.
% REF:	  Peterson & Abatzoglou (2014)
%
% IN:     gsi_raw - Daily GSI in row/column vector as class-type double.
%		  gsi_max - Maximum GSI value used to normalize daily GSI.
%		  N_DAYS - Number of days as a scalar value.
%		  N_YRS - Number of years as a scalar value.
%		  SCN - Scenario number between 1 and 3 indicating outer script 
%				progress; used to calculate or use gsi_max.
% OUT:    gu - Green-up as day of year with array dim [day,year].
%		  gsi_max - Max GSI to be used in normalizing daily raw GSI.
% CALL:   calcMovingAverage.m.
%==============================================================================


function [gu,gsi_max] = findGreenup(gsi_raw,gsi_max,N_DAYS,N_YRS,SCN);


%%=============================================================================
% Constants and initialized variables.
%==============================================================================
SEMI_DAY_WINDOW = 10;
gsi_mvg_avg = NaN(N_DAYS,N_YRS);


%%=============================================================================
% Moving average.
%==============================================================================
% Iterate over years then cast as single.
for yr = 1:N_YRS
	gsi_mvg_avg(:,yr) = calcMovingAverage(gsi_raw(:,yr),SEMI_DAY_WINDOW,N_DAYS);
end
gsi_mvg_avg = single(gsi_mvg_avg);


%%=============================================================================
% Normalize GSI.
%==============================================================================
% Find 95th percentile of the moving average for normalization.
if SCN ==1
	gsi_max = prctile(gsi_mvg_avg(:),95);
else
	gsi_max = gsi_max;
end

% Divide average GSI values by 95th percentile maximum to normalize.
gsi_norm = gsi_mvg_avg / gsi_max;


%%=============================================================================
% Green-up.
%==============================================================================
% Iterate over years to find green-up.
for yr=1:N_YRS

	% Find earliest day of year where GSI meets green-up threshold.
	greenup = min(find(gsi_norm(:,yr) >= 0.5));
	
	% Fill array, switching on empty.
    if ~isempty(greenup)
		gu(:,yr) = greenup;
    else
		gu(:,yr) = NaN;
    end 	% If statement.

end 		% i; N_YRS.



end 		% Function.