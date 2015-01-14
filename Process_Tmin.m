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

% Create model, experiment, and variable strings to be concatenated.
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

% Create path suffix and prefix strings to be concatenated.
FILE_PREFIX = '/storage/DOWNSCALED/CMIP5/MACAv2-METDATA/maca_v2_metdata_2var_10pat_CONUS_';
FILE_SUFFIX = {'_historical_tasmin.mat';...
               '_rcp45_tasmin.mat';...
               '_rcp85_tasmin.mat'};

% Break CONUS grid into regional subsets. The more powerful the machine, the
% larger the regional subsets can be.
LAT_START = [1];
LAT_END = [585];
LON_START = [1 701];
LON_END = [700 1386];

% Create constant for number of years, lat, lon, and models.
N_LAT = 585;
N_LON = 1386;
N_MDL = 20;
N_EXP = 3;

% Create and preallocate file.
out_file = matfile('MACA_Tmin.mat','writable',true);
out_file.data = NaN(N_LAT,N_LON,N_EXP,N_MDL,'single');

% main loop.
disp('Entering model loop...')
for i=1:N_MDL

	disp('Entering experiment loop...')
	for j=1:N_EXP

		% Concatenate prefixes and suffixes to create file name.
		file_name = [FILE_PREFIX,char(MDL_NAME(i)),char(FILE_SUFFIX(j))];
		file = matfile(file_name);
		disp(['Accessed ',char(MDL_NAME(i)),char(FILE_SUFFIX(j))])

		% Create n_yrs variable based on experiment.
		if j == 1
			n_yrs = 56;
		else
			n_yrs = 30;
		end

		% Load tmin data for days 60-151 (March-May).
		disp('Entering longitude loop and creating subsets...')
		for k=1:length(LON_START)

			% Break CONUS into regional lon subset.
            lon_subset = [LON_START(k):LON_END(k)];
            n_lon = length(lon_subset);

            % Preallocate model mean temperatures.
            disp('Preallocating variables...')
        	tmin_daily = NaN(n_yrs,N_LAT,n_lon,'single');
       		tmin_climo = NaN(N_LAT,n_lon,'single');

            % Pull years based on j (experiment). If historic, use all
            % years; otherwise use 2040 to 2069.
            disp('Loading TMIN data for days 60-151...')
            if j == 1
            	t_tmin = file.data(60:151,1:56,:,lon_subset);
            else
            	t_tmin = file.data(60:151,35:64,:,lon_subset);
            end

            % Calculate mean tmin for day window.
            disp('Taking daily mean for all years...')
            tmin_daily = squeeze(nanmean(t_tmin,1));

            % Calculate mean tmin over all years.
            disp('Taking climatological mean...')
            tmin_climo = squeeze(nanmean(tmin_daily,1));

            % Write output to file.
            disp('Writing output to file...')
            out_file.data(:,lon_subset,j,i) = tmin_climo;

        end 	% k; 1:LON_START
    end 	% j; 1:N_EXP
end 	% i; 1:N_MDL

% Exit MATLAB.
disp('Exiting MATLAB...')
exit