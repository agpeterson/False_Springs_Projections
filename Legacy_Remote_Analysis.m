%%=============================================================================
% NAME:   Model_Comparisons.m
% AUTHOR: Alexander Peterson
% DATE:   1 Dec. 2014
%
% DESC:   
% REF:	  None.
% NOTE:	  
%
% IN:     
% OUT:    
% CALL:   
%==============================================================================


%%=============================================================================
% Access LSF and GSI model output.
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

% Create years variables.
YEARS = [1950:2005 2006:2099 2006:2099];

% Years to average over.
hst_start = 1950;
hst_end = 2005;

fut1_start = 2040;
fut1_end = 2069;

fut2_start = 2070;
fut2_end = 2099;

% Find start and end indices.
hst_start_ind = find(YEARS == hst_start);
hst_end_ind = find(YEARS == hst_end);

fut1_start_ind = find(YEARS == fut1_start);
fut1_end_ind = find(YEARS == fut1_end);

fut2_start_ind = find(YEARS == fut2_start);
fut2_end_ind = find(YEARS == fut2_end);

% Create RCP/year indices.
hst_ind = hst_start_ind:hst_end_ind;

r45_fut1 = fut1_start_ind(1):fut1_end_ind(1);
r45_fut2 = fut2_start_ind(1):fut2_end_ind(1);

r85_fut1 = fut1_start_ind(2):fut1_end_ind(2);
r85_fut2 = fut2_start_ind(2):fut2_end_ind(2);

% Create structure for climatology indices.
clm_ind = {hst_ind r45_fut1 r45_fut2 r85_fut1 r85_fut2};

% Create constants for number of years, lat, lon, variables, and models.
N_LAT = 585;
N_LON = 1386;
N_MDL = 20;
N_VAR = 2;
N_YRS = 244;
N_CLM = 5;		% Hst + R45 fut1 + R45 fut2 + R85 fut1 + R85 fut2


%%=============================================================================
% Main loop to access file and derive metrics.
%==============================================================================

% Preallocate GU, LSF, and FSEI matfiles.
disp('Preallocating FS file...')
fsei_file = matfile('FSEI.mat','Writable',true);
fsei_file.clm_mean = NaN(N_CLM,N_LAT,N_LON,N_MDL,'single');

% Iterate over models for both variables.
for i=1:N_MDL

	% Concatenate to create filename for pointer.
	gu_file = [PATH_PREFIX,char(VARIABLE{1}),char(MODEL{i}),FILE_EXT];
	lsf_file = [PATH_PREFIX,char(VARIABLE{2}),char(MODEL{i}),FILE_EXT];
	
	% Create file pointers.
	gu_data = matfile(gu_file);
	lsf_data = matfile(lsf_file);

	% Print file to screen.
	disp(['Accessed file ', gu_file])
	disp(['Accessed file ', lsf_file])

	% Load model data into memory.
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

	% Preallocate FSEI variables.
	disp('Preallocating FSEI variables...')
	fs_sum = NaN(N_CLM,N_LAT,N_LON);
	gu_sum = NaN(N_CLM,N_LAT,N_LON);
	fsei = NaN(N_CLM,N_LAT,N_LON);

	% Derive FSEI.
	disp('Deriving FSEI...')
	for m=1:N_CLM
		fs_sum(m,:,:) = nansum(fs(clm_ind{m},:,:),1);	% Sums ones.
		gu_sum(m,:,:) = sum(~isnan(gu(clm_ind{m},:,:)),1);	% Sums # of elements.
		fsei(m,:,:) = (fs_sum(m,:,:) ./ gu_sum(m,:,:)) * 100;
	end

	% Write output to file.
	disp('Writing output to file...')
	fsei_file.clm_mean(:,:,:,i) = single(fsei);

	% Clear variables.
	clear gu_file lsf_file gu_data lsf_data fs fs_sum lsf_sum fsei
	clear lsf gu

end 	% j; 1:N_MDL


%%=============================================================================
% Derive climatological means and differences.
%==============================================================================

disp('Preallocating GU file...')
gu_file = matfile('GU.mat','Writable',true);
gu_file.clm_mean = NaN(N_CLM,N_LAT,N_LON,N_MDL,'single');
gu_file.clm_diff = NaN(N_CLM,N_LAT,N_LON,N_MDL,'single');

disp('Preallocating LSF file...')
lsf_file = matfile('LSF.mat','Writable',true);
lsf_file.clm_mean = NaN(N_CLM,N_LAT,N_LON,N_MDL,'single');
lsf_file.clm_diff = NaN(N_CLM,N_LAT,N_LON,N_MDL,'single');

% Iterate over models for both variables.
for i=1:N_MDL

	% Concatenate to create filename for pointer.
	gu_file = [PATH_PREFIX,char(VARIABLE{1}),char(MODEL{i}),FILE_EXT];
	lsf_file = [PATH_PREFIX,char(VARIABLE{2}),char(MODEL{i}),FILE_EXT];
	
	% Create file pointers.
	gu_data = matfile(gu_file);
	lsf_data = matfile(lsf_file);

	% Print file to screen.
	disp(['Accessed file ', gu_file])
	disp(['Accessed file ', lsf_file])

	% Load model data into memory.
	gu = double(gu_data.gsi_CONUS);
	lsf = double(lsf_data.lsf_CONUS);


	%%=========================================================================
	% Calculate climatological means and differences.
	%==========================================================================

	% Preallocate GU, LSF, and FSEI variables.
	disp('Preallocating GU and LSF variables...')
	gu_clm_mean = NaN(N_CLM,N_LAT,N_LON);
	gu_clm_diff = NaN(N_CLM,N_LAT,N_LON);
	lsf_clm_mean = NaN(N_CLM,N_LAT,N_LON);
	lsf_clm_diff = NaN(N_CLM,N_LAT,N_LON);

	% Iterate over N_CLM to aggregate.
	disp('Calculating climatological means and differences...')
	for j=1:N_CLM

		% Calculate climatological mean.
		gu_clm_mean(j,:,:) = nanmean(gu(clm_ind{j},:,:),1);
		lsf_clm_mean(j,:,:) = nanmean(lsf(clm_ind{j},:,:),1);

		% Take difference between historical and future periods.
		gu_clm_diff(j,:,:) = gu_clm_mean(j,:,:) - gu_clm_mean(1,:,:);
		lsf_clm_diff(j,:,:) = lsf_clm_mean(j,:,:) - lsf_clm_mean(1,:,:);
		
	end

	% Write output to file.
	disp('Writing output to files...')
	gu_file.clm_mean(:,:,:,i) = single(gu_clm_mean);
	gu_file.clm_diff(:,:,:,i) = single(gu_clm_diff);

	lsf_file.clm_mean(:,:,:,i) = single(lsf_clm_mean);
	lsf_file.clm_diff(:,:,:,i) = single(lsf_clm_diff);

	% Clear variables.
	clear gu_file lsf_file gu_data lsf_data lsf gu

end 	% j; 1:N_MDL


% Do the same as above but for FSEI differences.
disp('Preallocating FSEI file...')
fsei_file = matfile('FSEI.mat','Writable',true);
fsei_file.clm_diff = NaN(N_CLM,N_LAT,N_LON,N_MDL,'single');
fsei_file.clm_chng = NaN(N_CLM,N_LAT,N_LON,N_MDL,'single');

% Iterate over models.
for i=1:N_MDL

	% Load FSEI data into memory.
	disp('Loading model FSEI data into memory...')
	fsei = double(fsei_file.clm_mean(:,:,:,i));
	fsei_hst = squeeze(fsei(1,:,:));


	% Preallocate per model.
	disp('Preallocating variables per model...')
	fsei_diff = NaN(N_CLM,N_LAT,N_LON);
	fsei_chng = NaN(N_CLM,N_LAT,N_LON);

	% Iterate over N_CLM to aggregate.
	disp('Calculating FSEI differences and changes...')
	for j=1:N_CLM

		% Take difference between historical and future periods.
		fsei_diff(j,:,:) = squeeze(fsei(j,:,:)) - fsei_hst;

		% Calculate percent change between periods.
		fsei_chng(j,:,:) = ((squeeze(fsei(j,:,:)) - fsei_hst) ./ ...
			fsei_hst) * 100;

	end 	% j; 1:N_CLM

	% Write output to file.
	disp('Writing output to file...')
	fsei_file.clm_diff(:,:,:,i) = single(fsei_diff);
	fsei_file.clm_chng(:,:,:,i) = single(fsei_chng);

end 	% i; 1:N_MDL


%%=============================================================================
% Explore driving changes in GU.
%==============================================================================

% Load data into memory.
load('MACA_Tmin.mat');

% Convert from K to C.
data = data - 273.15;

% Create variables for historic and future periods.
tmin_hst = squeeze(data(:,:,1,:));
tmin_fut = squeeze(data(:,:,3,:));

% Take difference for delta.
tmin_delta = tmin_fut - tmin_hst;

% Create multi-model mean delta C.
tmin_mdl_delta = squeeze(nanmean(tmin_delta,3));

% Save delta to be used in sensitivity experiment.
save('tmin_delta.mat','tmin_mdl_delta')


%%=============================================================================
% Sensitivity analysis/experiments.
%==============================================================================

gu_exp_file = matfile('GU_Sensitivity.mat');
lsf_exp_file = matfile('LSF_Sensitivity.mat');

% Set up variables and data before calling calcFSEI() function.
gu = gu_exp_file.data;
lsf = lsf_exp_file.data;
day_lag = 7;
yrs_ind = 1:30;
climo_ind = {yrs_ind}

% Derive FSEI.
[fs_exp,fsei_exp] = calcFSEI(gu,lsf,day_lag,climo_ind);

% Convert to single and save.
fsei = single(fsei_exp);           
fs = single(fs_exp);
save('FSEI_Sensitivity.mat','fsei','fs')

% Determine significance between experiments and observed.
gu_obs_file = matfile('GU_METDATA')
lsf_obs_file = matfile('LSF_METDATA')

% Load data into memory.
gu_obs = gu_obs_file.data;
lsf_obs = lsf_obs_file.data;

% Derive FSEI.
[fs_obs,fsei_obs] = calcFSEI(gu_obs,lsf_obs,day_lag,climo_ind);

% Convert to single and save.
fsei = single(fsei_obs);
fs = single(fs_obs);
save('FSEI_METDATA.mat','fsei','fs')

% Determine change and statistical significance for each model.
% Create constants for lat, lon, and models.
N_LAT = 585;
N_LON = 1386;
N_MDL = 20;
N_YRS = 30;

% Open matfile and preallocate.
file = matfile('FSEI_METDATA_Change.mat','Writable',true);
file.fs_change = NaN(N_LAT,N_LON,N_MDL,'single');
file.fs_significance = NaN(N_LAT,N_LON,N_MDL,'single');

% Year indices variable.
yr_ind = 1:30;

% Set parallel pool options.
opts = statset('UseParallel','always');

% Open parallel processor pool.
disp('Opening parallel processor pool...')
if matlabpool('size') == 0
    matlabpool open local 12
end

% Iterate over all sensitivity runs.
for i=1:N_MDL
	tic
	disp({'Model: ' i})

	% Create temporary variable for each sensitivity run.
	disp('Creating temporary variable for each sensitivity run...')
	fs_tmp = squeeze(fs_exp(:,:,:,i));

	% Bootstrap resampling. Iterate over lon and call bootstrp() function using
	% 1000 samples, taking the mean of the difference between the future and 
	% historic periods.
	disp('Preallocating bootstrap variable...')
	fs_bootstrap = NaN(1000,N_LAT,N_LON);

	disp('Bootstrap function call...')
	for x=1:N_LON
	    fs_bootstrap(:,:,x) = bootstrp(1000,@nanmean,...
	    	fs_tmp(yr_ind,:,x)-fs_obs(yr_ind,:,x),'Options',opts);
	end

	% Change and significance.
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

	% Write output to file.
	disp('Writing output to file...')
	file.fs_change(:,:,i) = single(fs_change);
	file.fs_significance(:,:,i) = single(fs_significance);

	toc
	% Clear variables.
	clear fs_tmp fs_bootstrap fs_greater fs_less fs_increase fs_decrease
	clear fs_significance fs_change

end 	% i; 1:N_MDL


%% Create multi-model mean FSEI change and significance.
% Take the mean of change, find where 90% of models agree on significance.

fs_mdl_change = squeeze(nanmean(double(fs_change),3));
test = nansum(fs_significance,3);

for x=1:1386
	for y=1:585
		if test(y,x) >= (20*.75)
			test2(y,x) = 1;
		elseif test(y,x) < (20*.75)
			test2(y,x) = 0;
		end
	end
end





