%%=============================================================================
%
% NAME:   FSEI_Change_Significance.m
% AUTHOR: Alexander Peterson
% DATE:   19 Jan. 2015
% DESC:   Determine statistical significance in false spring changes between
%			future and historic periods using a bootstrap resampling method.
%
%==============================================================================


%%=============================================================================
% Variables and constants.
%==============================================================================

% Path suffix and prefix strings to be concatenated for model access.
PATH_PREFIX = 'False_Springs_Data/';	% Thunder
VARIABLE = {'gu_';...
			'lsf_'};
MODEL = {'bcc-csm1-1';...
         'bcc-csm1-1-m';...
	     'BNU-ESM';...
	     'CanESM2';...
	     'CCSM4';...
	     'CNRM-CM5';...
	     'CSIRO-Mk3-6-0';...
	     'GFDL-ESM2G';...
	     'GFDL-ESM2M';...
	     'HadGEM2-CC365';...
	     'HadGEM2-ES365';...
	     'inmcm4';...
	     'IPSL-CM5A-LR';...
	     'IPSL-CM5A-MR';...
	     'IPSL-CM5B-LR';...
	     'MIROC5';...
	     'MIROC-ESM-CHEM';...
	     'MIROC-ESM';...
	     'MRI-CGCM3';...
	     'NorESM1-M'};
FILE_EXT = '.mat';

% Create constants for lat, lon, and models.
N_LAT = 585;
N_LON = 1386;
N_MDL = 20;
N_YRS = 244;

% Create indices for historic and future periods. Use a random set of integers
% as historic indices.
hst_ind = randi([1 56],1,30);
fut_ind = 185:214;

% Set parallel pool options.
opts = statset('UseParallel','always');

% Open parallel processor pool.
disp('Opening parallel processor pool...')
if matlabpool('size') == 0
    matlabpool open local 12
end


%%=============================================================================
% Main loop - access LSF and GU model output to calculate false springs, then
% determine statistical significance using bootstrap resampling.
%==============================================================================

% Open matfile and preallocate.
file = matfile('FSEI_Change_Significance.mat','Writable',true);
file.fs_change = NaN(N_LAT,N_LON,N_MDL,'single');
file.fs_significance = NaN(N_LAT,N_LON,N_MDL,'single');

% Iterate over models for both variables.
for i=1:N_MDL

	% Concatenate to create filename for pointer.
	gu_file = [PATH_PREFIX,char(VARIABLE{1}),char(MODEL{i}),FILE_EXT];
	lsf_file = [PATH_PREFIX,char(VARIABLE{2}),char(MODEL{i}),FILE_EXT];
	
	% Create file pointers.
	gu_data = matfile(gu_file);
	lsf_data = matfile(lsf_file);

	% Print file to screen.
	disp(['Accessed model ', gu_file])
	disp(['Accessed model ', lsf_file])

	% Load model data into memory.
	disp('Loading data into memory...')
	gu = double(gu_data.gsi_CONUS);
	lsf = double(lsf_data.lsf_CONUS);


	%%=========================================================================
	% Derive false springs.
	%==========================================================================

	% Preallocate FS variable.
	disp('Preallocating FS variable...')
	fs = NaN(N_YRS,N_LAT,N_LON);

	% Iterate over lon, lat, and years.
	disp('Classifying false springs...')
	for j=1:N_LON
		for k=1:N_LAT
			for l=1:N_YRS
				if isnan(lsf(l,k,j)) || isnan(gu(l,k,j))
					fs(l,k,j) = NaN;
				elseif lsf(l,k,j) >= (7 + gu(l,k,j))
					fs(l,k,j) = 1;
				else
					fs(l,k,j) = 0;
				end 	% If statement.
			end 	% l; 1:N_YRS
		end 	% k; 1:N_LAT
	end 	% j; 1:N_LON


	%%=========================================================================
	% Bootstrap resampling.
	%==========================================================================

	% Iterate over lon and call bootstrp() function using 1000 samples, taking 
	% the mean of the difference between the future and historic periods.
	disp('Preallocating bootstrap variable...')
	fs_bootstrap = NaN(1000,N_LAT,N_LON);

	disp('Bootstrap function call...')
	for x=1:N_LON
	    fs_bootstrap(:,:,x) = bootstrp(1000,@nanmean,...
	    	fs(fut_ind,:,x)-fs(hst_ind,:,x),'Options',opts);
	end


	%%=========================================================================
	% Change and significance.
	%==========================================================================

	% Magnitude of change.
	disp('Calculating magnitude of change...')
	fs_change = squeeze(nanmean(fs_bootstrap,1));

	% Find lat/lon greater/less than 0 to determine increasing/decreasing
	% false springs.
	disp('Determining statistical significance...')
	fs_greater = squeeze(nansum(fs_bootstrap > 0));
	fs_less = squeeze(nansum(fs_bootstrap < 0));
	fs_increase = fs_greater >= 0.95*1000;
	fs_decrease = fs_less >= 0.95*1000;

	% Iterate over lat/lon and classify significance.
	for x=1:N_LON
	    for y=1:N_LAT
	        if fs_decrease(y,x) == 1 || fs_increase(y,x) == 1
	            fs_significance(y,x) = 1;
	        else
	            fs_significance(y,x) = 0;
	        end
	    end
	end


	%%=========================================================================
	% Write output to file; close MATLAB when complete.
	%==========================================================================

	% Write output to file.
	disp('Writing output to file...')
	file.fs_change(:,:,i) = single(fs_change);
	file.fs_significance(:,:,i) = single(fs_significance);

	% Clear variables.
	clear gu lsf fs*

end 	% i; 1:N_MDL

% Close parallel pool and program.
disp('Script complete; closing parallel pool and MATLAB...')
matlabpool close;
exit