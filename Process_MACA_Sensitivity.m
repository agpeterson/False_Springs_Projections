%%============================================================================
%
% NAME:   Process_MACA_Sensitivity.m
% AUTHOR: Alexander Peterson
% DATE:   27 Jan. 2015
% DESC:   Script containing code analyzing GU, LSF, and FSEI.
%
%=============================================================================



%%============================================================================
% Initial variables and constants.
%-----------------------------------------------------------------------------

% Create path suffix and prefix strings to be concatenated.
PATH = '/storage/DOWNSCALED/CMIP5/MACAv2-METDATA/';
FILE_PREFIX = 'maca_v2_metdata_2var_10pat_CONUS_';
FILE_SUFFIX = '_historical_tasmin.mat';

% Model names.
MDL_NAME = {'bcc-csm1-1';...
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


%-----------------------------------------------------------------------------
% Create constants for number of years, lat, lon, variables, and models.
%-----------------------------------------------------------------------------

% Break CONUS grid into regional subsets. The more powerful the machine, the
% larger the regional subsets can be.
LAT_START = [1];
LAT_END = [585];
LON_START = [1 701];
LON_END = [700 1386];

% Create constant for number of years, lat, lon, and models.
N_LAT = 585;
N_LON = 1386;
N_YRS = 56;
N_DAYS = 181;
N_MDL = 20;

% Create time indices.
DAY_IND = 1:181;

% Load tmin deltas for sensitivity experiment.
disp('Loading MACA tmin delta...')
load('Tmin_Delta.mat');



%%============================================================================
% GSI constants.
%-----------------------------------------------------------------------------

disp('Creating day length and VPD variables...')

% Open day_length data.
load('Day_Length_MACA.mat')
day_length = day_length(DAY_IND,:);

% Find where photoperiod is greater than or less than bounds and normalize; 
% all values below bounds set to 0 and all values above bounds set to 1.
HOUR_LOW = 10;
HOUR_HIGH = 11;
day_length = day_length(DAY_IND,:);
dayl = NaN(size(day_length)); 
f = find(day_length > HOUR_LOW & day_length < HOUR_HIGH);
dayl(f) = day_length(f) - HOUR_LOW;
dayl(day_length <= HOUR_LOW) = 0;
dayl(day_length >= HOUR_HIGH) = 1;

% VPD. Assume held constant across time and space.
VPD_LOW = 0.9;
VPD_HIGH = 4.1;
vpd = single(1.0 - (1.0 - VPD_LOW) / (VPD_HIGH-VPD_LOW));



%%============================================================================
% Create modeled climatological mean GU, LSF, and FSEI.
%-----------------------------------------------------------------------------

% Open parallel processor pool.
if matlabpool('size') == 0
    matlabpool open local 12
end

% Create matfile and preallocate GU, LSF, and FSEI values.
disp('Preallocating Model_Sensitivity.mat file...')
mdl_senstvty_file = matfile('Model_Sensitivity.mat','Writable',true);
mdl_senstvty_file.gu = NaN(N_YRS,N_LAT,N_LON,N_MDL,'single');
mdl_senstvty_file.lsf = NaN(N_YRS,N_LAT,N_LON,N_MDL,'single');
mdl_senstvty_file.fs = NaN(N_YRS,N_LAT,N_LON,N_MDL,'single');
mdl_senstvty_file.gu_mean = NaN(N_LAT,N_LON,N_MDL,'single');
mdl_senstvty_file.lsf_mean = NaN(N_LAT,N_LON,N_MDL,'single')
mdl_senstvty_file.fsei = NaN(N_LAT,N_LON,N_MDL,'single');

% Iterate over models for both variables, derive false springs, and take
% climatological means.
disp('Iterating over models to derive false springs and climatologies...')
for mdl=1:N_MDL

	% Concatenate to create filename for pointer.
	disp(['Accessing historical tmin data from ' char(MDL_NAME{mdl}) '...'])
	tmin_file = matfile([PATH,FILE_PREFIX,char(MDL_NAME{mdl}),FILE_SUFFIX]);

	% Iterate through regions.
	for x=1:length(LON_START)
	        
	    % Break CONUS into regional lon subset.
	    lon_subset = [LON_START(x):LON_END(x)];
	    n_lon = length(lon_subset);

	    for y=1:length(LAT_START)

	        % Write lon and lat cell to output.
	        disp({'Lon: ',x; 'Lat: ',y})
	    
	        % Break CONUS into regional lat subsets.
	        lat_subset = [LAT_START(y):LAT_END(y)];
	        n_lat = length(lat_subset);

	        % Preallocate subset variables.
	        gu_subset = NaN(N_YRS,n_lat,n_lon,'single');
	        lsf_subset = NaN(N_YRS,n_lat,n_lon,'single');

	        % Break day length into regional subset.
	        dayl_subset = dayl(:,lat_subset);

	        % Create temporary variable to store daily and yearly data for
	        % each lat/lon subset. Replace -9999 (missing values) with NaN.
	        disp('Loading temperature data...')
	        t_tmin = tmin_file.data(DAY_IND,:,lat_subset,lon_subset);

	        % Extract tmin deltas for lat/lon subsets.
            t_delta = squeeze(tmin_delta(lat_subset,lon_subset,mdl));
            t_delta = flipud(t_delta);

            % Parallel iteration over lon_subset.
            tic
            disp('Entering parallel loop...')
            parfor lon=1:n_lon

                % Call parallel function.
                [lsf_sub,gu_sub] = findSpringEvents(t_tmin(:,:,:,lon),...
                                                    t_delta(:,lon),...
                                                	dayl_subset(:,:),...
                                                    vpd,N_DAYS,N_YRS,n_lat,1);

                % Concatenate subregional variables to subset.
                gu_subset(:,:,lon) = gu_sub;
                lsf_subset(:,:,lon) = lsf_sub;

            end     % parfor - lon; n_lon.
            toc

            % Classify false springs and derive FSEI.
			disp('Calling calcFSEI() function...')
			[fs,fsei] = calcFSEI(gu_subset,lsf_subset,7);

			% Calculate climatological mean GU and LSF.
			disp('Taking climatological mean GU and LSF...')
			gu_mean = squeeze(nanmean(gu_subset,1));
			lsf_mean = squeeze(nanmean(lsf_subset,1));

			% Write regional subsets to file.
            disp('Writing output to file...')
            mdl_senstvty_file.gu(:,lat_subset,lon_subset,mdl) = single(lsf_subset);
            mdl_senstvty_file.lsf(:,lat_subset,lon_subset,mdl) = single(gu_subset);
            mdl_senstvty_file.fs(:,lat_subset,lon_subset,mdl) = single(fs);
			mdl_senstvty_file.gu_mean(lat_subset,lon_subset,mdl) = single(gu_mean);
			mdl_senstvty_file.lsf_mean(lat_subset,lon_subset,mdl) = single(lsf_mean);
			mdl_senstvty_file.fsei(lat_subset,lon_subset,mdl) = single(fsei);

			% Clear variables.
			clear gu_subset lsf_subset dayl_subset t_tmin t_delta lsf_sub gu_subset
			clear fs fsei gu_mean lsf_mean

		end 	% y; 1:length(LAT_START)
	end 		% x; 1:length(LON_START)
end 			% mdl; 1:length(N_MDL)

% Exit MATLAB.
exit