%==============================================================================
% NAME:   calculategsi.m
% DESC:   As per Jolly et al., 2005, this function calculates GSI using minimum
%		  temperature, photoperiod (day length), and VPD as proxy for
%		  evaporative demand. Each component is run through a step function,
%		  where 0 is constrained and 1 is optimal.
% IN:     Daily minimum temperature, photoperiod, and VPD data.
% OUT:    GSI.
% CALL:   None.
% AUTH:   John Abatzoglou. Modified by Alexander Peterson, 21 Sept. 2014.
% NOTE:	  Needs further refactoring and clean-up.
%==============================================================================

function [gsi] = calcgsi(tmin,photoperiod,vpd)

% Photoperiod. Find where photoperiod is greater than or less than bounds and
% normalize; all values below bounds set to 0 and all values above bounds set
% to 1.
HOUR_LOW = 10;
HOUR_HIGH = 11;
f = find(photoperiod > HOUR_LOW & photoperiod < HOUR_HIGH);
dayl(f) = (photoperiod(f)-HR_LOW) /HR_HIGH ;
dayl(photoperiod <= HR_LOW) = 0;
dayl(photoperiod >= HR_HIGH) = 1;
dayl = dayl';

% Minimum temperature. Follow above.
TMIN_LOW = -2;
TMIN_HIGH = 5;
f1 = find(tmin > T_LOW & tmin < T_HIGH);
tmin(f1) = (tmin(f1)-T_LOW) /T_HIGH;
tmin(tmin <= T_LOW) = 0;
tmin(tmin >= T_HIGH) = 1;

% VPD. For now, we assume VPD is one for all lat/lon. Before actual use, this
% needs to be double-checked for errors, especially regarding lower and upper
% bounds.
% vpd(vpd <= .9) = 1;
% f = find(vpd > .9 & vpd < 4.1);
% vpd(f) = 1-(vpd(f)-.9) / 3.2;
% vpd(vpd >= 4.1) = 0;
VPD_LOW = 0.9;
VPD_HIGH = 4.1;
f1 = find(vpd > V_LOW & vpd < V_HIGH);
vpd(f1) = 1.0 -(vpd(f1)-V_LOW) / V_HIGH;
vpd(vpd <= V_LOW) = 1;
vpd(vpd >= V_HIGH) = 0;

% Reshape dayl and calculate GSI.
dayl = reshape(dayl,size(vpd,1),size(vpd,2));
gsi = tmin.*vpd.*dayl;

end 	% Function.