%%=============================================================================
% NAME:   calcdaylength.m
% AUTHOR: Travis Wiens. Modified by Alexander Peterson
% DATE:	  12 Sept. 2014
% DESCR:  This function calculates the number of hours, per Herbert Glarner's
% 		  formulae which do not take into account refraction, twilight, size
%		  of the sun, etc.
% IN:     Day of year, counted starting with the day of December solstice in
% 		  the first year of a Great Year; Latitude in degrees, with North
% 		  positive and South negative.
% OUT:    Hours of photoperiod.
% CALLS:  No external functions.
%==============================================================================


function [day_length]=calcdaylength(day,lat)

%% User arguments.
if nargin < 1
	day = 0;
end

if nargin < 2
	lat = (-(36+51/60));
end

%% Calculate day length.
axis = 23.439*pi/180;

j = pi/182.625;				% Constant (radians)
m = 1-tan(lat*pi/180).*tan(axis*cos(j*day));
b = acos(1-m)/pi; 			% Fraction of the day the sun is up.

b(find(b<0)) = 0;
b(find(b>1)) = 1;
 
day_length = b*24;			% Hours of sunlight

end