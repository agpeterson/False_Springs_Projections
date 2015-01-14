%%=============================================================================
% NAME:   calcGSI.m
% AUTHOR: John Abatzoglou, Alexander Peterson, Katherine Hegewisch
% DATE:   1 Nov. 2014
%
% DESC:   This function calculates the Growing Season Index (GSI) using minimum 
%		  temperature, photoperiod (day length), and vapor pressure deficit 
%		  (VPD) as proxy for evaporative demand. Each component is run through 
%		  a step function, where 0 is constrained and 1 is optimal.
% REF:	  Jolly et al. (2005)
% NOTE:	  Inputs can be class-type single or double.
%
% IN:     tmin - Daily minimum temperature in units C in row/column vector.
%		  photoperiod - Photoperiod for all calendar days in row/column vector.
%		  vpd - Vapor pressure deficit in row/column vector.
% OUT:    gsi - GSI in row/column vector.
% CALL:   None.
%==============================================================================


function [gsi] = calcGSI(tmin,photoperiod,vpd)


%%=============================================================================
% Squeeze inputs.
%==============================================================================
tmin = squeeze(tmin);
photoperiod = squeeze(photoperiod);
vpd = squeeze(vpd);


%%=============================================================================
% Minimum temperature.
%==============================================================================
TMIN_LOW = -2;
TMIN_HIGH = 5;
temp = NaN(size(tmin),'single');
f1 = find(tmin > TMIN_LOW & tmin < TMIN_HIGH);
temp(f1) = (tmin(f1)-TMIN_LOW) /(TMIN_HIGH-TMIN_LOW);
temp(tmin <= TMIN_LOW) = 0;
temp(tmin >= TMIN_HIGH) = 1;


%%=============================================================================
% Photoperiod.
% *NOTE: For now, photoperiod iis calculated outside and passed into the
% function to reduce processing time.
%==============================================================================
% HOUR_LOW = 10;
% HOUR_HIGH = 11;
% f = find(photoperiod > HOUR_LOW & photoperiod < HOUR_HIGH);
% dayl(f) = (photoperiod(f)-HOUR_LOW) /HOUR_HIGH ;
% dayl(photoperiod <= HOUR_LOW) = 0;
% dayl(photoperiod >= HOUR_HIGH) = 1;


%%=============================================================================
% VPD.
%==============================================================================
% VPD_LOW = 0.9;
% VPD_HIGH = 4.1;
% f1 = find(vpd > VPD_LOW & vpd < VPD_HIGH);
% vpd(f1) = 1.0 -(vpd(f1)-VPD_LOW) / VPD_HIGH;
% vpd(vpd <= VPD_LOW) = 1;
% vpd(vpd >= VPD_HIGH) = 0;


%%=============================================================================
% GSI.
%==============================================================================
gsi = vpd .* temp .* photoperiod;



end 	% Function.