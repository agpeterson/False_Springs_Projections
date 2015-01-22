%%=============================================================================
%
% NAME:   calcFSEI.m
% AUTHOR: Alexander Peterson, John Abatzoglou
% DATE:   19 Jan. 2015
% DESC:   
% REF:	  Peterson & Abatzoglou (2014)
% IN:     
% OUT:    
% CALL:   
%
%==============================================================================

function [fs,fsei] = calcFSEI(gu,lsf,DAY_LAG,CLIMO_IND)


%%=============================================================================
% Constants and initialized variables.
%==============================================================================

N_YRS = size(gu,1);
N_LAT = size(gu,2);
N_LON = size(gu,3);
N_CLM = size(CLIMO_IND);

% If GU has 4 dimensions, use size as N_MDL.
if ndims(gu) == 4
	N_MDL = size(gu,4);
else
	N_MDL = 0;
end


%%=============================================================================
% Main loop - switch on N_MDL to determine variable sizes.
%==============================================================================

disp('Entering if statement...')
if N_MDL == 0

	% Preallocate FS variable.
	disp('Preallocating FS variable...')
	fs = NaN(N_YRS,N_LAT,N_LON);

	% Iterate over lon, lat, and years.
	disp('Classifying false springs...')
	for lon=1:N_LON
		for lat=1:N_LAT
			for yrs=1:N_YRS
				if isnan(lsf(yrs,lat,lon)) || isnan(gu(yrs,lat,lon))
					fs(yrs,lat,lon) = NaN;
				elseif lsf(yrs,lat,lon) >= (DAY_LAG + gu(yrs,lat,lon))
					fs(yrs,lat,lon) = 1;
				else
					fs(yrs,lat,lon) = 0;
				end 	% if statement
			end 	% yrs; 1:N_YRS
		end 	% lat; 1:N_LAT
	end 	% lon; 1:N_LON

	% Check for number of climatologies to create.
	if N_CLM == 1

		% Preallocate FSEI variables.
		disp('Preallocating FSEI variables...')
		fs_sum = NaN(N_LAT,N_LON);
		gu_sum = NaN(N_LAT,N_LON);
		fsei = NaN(N_LAT,N_LON);

		% Derive FSEI.
		disp('Deriving FSEI...')
		fs_sum(:,:) = nansum(fs(CLIMO_IND{1},:,:),1);	% Sums ones.
		gu_sum(:,:) = sum(~isnan(gu(CLIMO_IND{1},:,:)),1); % Sums # of elements.
		fsei(:,:) = (fs_sum(:,:) ./ gu_sum(:,:)) * 100;

	elseif N_CLM > 1

		% Preallocate FSEI variables.
		disp('Preallocating FSEI variables...')
		fs_sum = NaN(N_LAT,N_LON,N_CLM);
		gu_sum = NaN(N_LAT,N_LON,N_CLM);
		fsei = NaN(N_LAT,N_LON,N_CLM);

		% Derive FSEI.
		disp('Deriving FSEI...')
		for clm=1:N_CLM
			fs_sum(:,:,clm) = nansum(fs(CLIMO_IND{clm},:,:),1);	% Sums ones.
			gu_sum(:,:,clm) = sum(~isnan(gu(CLIMO_IND{clm},:,:)),1);	% Sums # of elements.
			fsei(:,:,clm) = (fs_sum(clm,:,:) ./ gu_sum(clm,:,:)) * 100;
		end

	end 	% if statement; N_CLM == 1 or N_CLM > 1

elseif N_MDL > 1

	% Preallocate FS variable.
	disp('Preallocating FS and FSEI variables...')
	fs = NaN(N_YRS,N_LAT,N_LON,N_MDL);
	
	% Switch on N_CLM.
	if N_CLM == 1
		fsei = NaN(N_LAT,N_LON,N_MDL);
	elseif N_CLM > 1
		fsei = NaN(N_LAT,N_LON,N_CLM,N_MDL);
	end

	% Iterate over models.
	for mdl=1:N_MDL

		disp({'Model: ' mdl})

		% Pull individual models.
		disp('Loading individual model data...')
		gu_tmp = squeeze(gu(:,:,:,mdl));
		lsf_tmp = squeeze(lsf(:,:,:,mdl));

		% Preallocate temporary fs variable.
		disp('Preallocating fs_tmp variable...')
		fs_tmp = NaN(N_YRS,N_LAT,N_LON);
		
		% Iterate over lat/lon/years to classify false springs.
		for lon=1:N_LON
			for lat=1:N_LAT
				for yrs=1:N_YRS
					if isnan(lsf_tmp(yrs,lat,lon)) || isnan(gu_tmp(yrs,lat,lon))
						fs_tmp(yrs,lat,lon) = NaN;
					elseif lsf(yrs,lat,lon) >= (DAY_LAG + gu_tmp(yrs,lat,lon))
						fs_tmp(yrs,lat,lon) = 1;
					else
						fs_tmp(yrs,lat,lon) = 0;
					end 	% if statement
				end 	% yrs; 1:N_YRS
			end 	% lat; 1:N_LAT
		end 	% lon; 1:N_LON
	
		% Check for number of climatologies.
		if N_CLM == 1

			disp('Preallocating FSEI variables...')
			fs_sum = NaN(N_LAT,N_LON);
			gu_sum = NaN(N_LAT,N_LON);
			fsei_tmp = NaN(N_LAT,N_LON);
			
			% Derive FSEI.
			disp('Deriving FSEI...')
			fs_sum = squeeze(nansum(fs_tmp(CLIMO_IND{1},:,:),1));	% Sums ones.
			gu_sum = squeeze(sum(~isnan(gu_tmp(CLIMO_IND{1},:,:)),1));	% Sums # of elements.
			fsei_tmp = (fs_sum(:,:) ./ gu_sum(:,:)) * 100;
			
			% Write fsei_tmp to fsei variable.
			disp('Writing *_tmp variables to FSEI and FS...')
			fs(:,:,:,mdl) = fs_tmp;
			fsei(:,:,mdl) = fsei_tmp;

			% Clear variables.
			clear *_sum *_tmp

		elseif N_CLM > 1

			disp('Preallocating FSEI variables...')
			fs_sum = NaN(N_LAT,N_LON,N_CLM);
			gu_sum = NaN(N_LAT,N_LON,N_CLM);
			fsei_tmp = NaN(N_LAT,N_LON,N_CLM);
			
			% Derive FSEI.
			disp('Deriving FSEI...')
			for clm=1:N_CLM
				fs_sum(:,:,clm) = nansum(fs(CLIMO_IND{clm},:,:),1);	% Sums ones.
				gu_sum(:,:,clm) = sum(~isnan(gu_tmp(CLIMO_IND{clm},:,:)),1);	% Sums # of elements.
				fsei_tmp(:,:,clm) = (fs_sum(clm,:,:) ./ gu_sum(clm,:,:)) * 100;
			end

			% Write fsei_tmp to fsei variable.
			disp('Writing *_tmp variables to FSEI and FS...')
			fs(:,:,:,mdl) = fs_tmp;
			fsei(:,:,:,mdl) = fsei_tmp;

			% Clear variables.
			clear *_sum *_tmp

		end 	% if statement; N_CLM == 1 or N_CLM > 1
	end 	% mdl; 1:N_MDL
end 	% if statement


end 	% Function.
