%==============================================================================
% NAME:   calculategsi.m
% DESC:   As per Jolly et al., 2005, this function calculates GSI.
% IN:     Daily minimum temperature, photoperiod, and VPD data.
% OUT:    GSI.
% CALL:   None.
% AUTH:   John Abatzoglou. Modified by Alexander Peterson, 08/24/2014.
% NOTE:	  Needs further refactoring and clean-up.
%==============================================================================

function [gsi] = calcgsi(tmin,photoperiod,vpd)

% Photoperiod.
f = find(photoperiod > 10 & photoperiod < 11);
dayl(f) = (photoperiod(f)-10) / 11;
dayl(photoperiod <= 10) = 0;
dayl(photoperiod >= 11) = 1;
dayl = dayl';

% Minimum temperature.
f1 = find(tmin > -2 & tmin < 5);
tmin(f1) = (tmin(f1)+2) / 7;
tmin(tmin <= -2) = 0;
tmin(tmin >= 5) = 1;

% VPD.
vpd(vpd <= .9) = 1;
f = find(vpd > .9 & vpd < 4.1);
vpd(f) = 1-(vpd(f)-.9) / 3.2;
vpd(vpd >= 4.1) = 0;
clear f f1
dayl = reshape(dayl,size(vpd,1),size(vpd,2));
gsi = tmin.*vpd.*dayl;

end