%%=============================================================================
% NAME:   Calculate_Day_Length.m
% AUTHOR: Alexander Peterson
% DATE:   16 Sept. 2014
% DESCR:  This script calls one model within the MACAv2-METDATA dataset and
%		  calculates day length (photoperiod) across all latitudes.
% IN:     MACAv2-METDATA
% OUT:    N/A
% CALLS:  calcdaylength.m
%==============================================================================

% Create file name string.
model = char(MDL_NAME(MDL_TARGET(1)));
file_name = [PATH_PREFIX,FILE_PREFIX,model,char(FILE_SUFFIX(1))];

% Create file pointer and pull latitude.
file = matfile(file_name);
lat = file.lat;

% Iterate over latitude to calculate day length for all days.
for i=1:N_LAT
	day_length(i,:) = calcdaylength(1:N_DAY,lat(i));
end
day_length = double(day_length');     % Transpose such that latidue is outside 
									  % dimension and change class to double.

% Save.
save('day_length.mat','day_length')