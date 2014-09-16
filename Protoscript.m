%%=============================================================================
% NAME:   Protoscript.m
% AUTHOR: Alexander Peterson
% DATE:   12 Sept. 2014
% DESCR:  This script contains prototype code to read through MACAv2-METDATA
%         products.
% IN:     MACAv2-METDATA
% OUT:    N/A
% CALLS:  findlsf.m
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
EXP_NAME = {'rcp45';...
            'rcp85'};
VAR_NAME = {'tasmax';...
            'tasmin';...
            'rhsmax';...
            'rhsmin';...
            'pr';...
            'rsds';...
            'uas';...
            'vas';...
            'huss'};

% Set target data by specifying index for above constants.
MDL_TARGET = [1:20];
EXP_TARGET = [1:2];
VAR_TARGET = [2];

% Create path suffix and prefix strings to be concatenated.
PATH_PREFIX = 'maca_v2_metdata_2var_10pat_CONUS_';
PATH_SUFFIX = {'_historical_tasmin.mat';...
               '_rcp45_tasmin.mat';...
               '_rcp85_tasmin.mat'};

% Break CONUS grid into regional subsets.
LAT_START = [1 200 400];
LAT_END = [199 399 585];
LON_START = [1 300 600 900 1100];
LON_END = [299 599 899 1099 1386];

% Create constant for number of years, lat, and lon.
N_YEARS = 151;
N_LAT = 585;
N_LON = 1386;

% Preallocate CONUS variables.
lsf_CONUS = NaN(N_YEARS,N_LAT,N_LON);
gsi_CONUS = NaN(N_YEARS,N_LON,N_LON);


%%=============================================================================
% Body of script to access and process MACA data. To do so, iterate over models
% and experiments to concatenate strings and access the .mat files, using
% nested loops to subset CONUS and process the data.
%==============================================================================

% Model and experiment iteration. 
for m=1:length(MDL_TARGET)
    
    % Subset model name using character array.
    model = char(MDL_NAME(MDL_TARGET(m)));

    % Concatenate component strings for file path over experiments.
    for e=1:3  % REFACTOR - should not be hard-coded.
		
        % Create path string for each file.
    	file_name = [PATH_PREFIX,model,char(PATH_SUFFIX(e))];
    		
    	% Set pointer as matfile.
    	m = matfile(file_name);

        % Create variable to hold number of years.
        n_years = size(m.data,2);

        % Extract lat values for one model and calculate day length.
        lat = m.lat;
        if m == 1
            for y=1:N_LAT
                day_length(y,:) = calcdaylength(1:366,m.lat(y));    % REFACTOR - check files before continuing with this code.
            end
        end

        %%=====================================================================
    	% Iterate over CONUS by breaking lat/lon into regional subset, then
        % iterating over regional subset to call findlsf function.
        %======================================================================
    	
        for x=1:length(LON_START)
        
            % Break CONUS into regional lon subset.
            lon_subset = [LON_START(x):LON_END(x)];

            for y=1:length(LAT_START)
                
                % Break CONUS into regional lat subset.
                lat_subset = [LAT_START(y):LAT_END(y)];
                day_subset = day_length(lat_subset,:);

                % Preallocate subset variables.
                lsf_subset = NaN(n_years,length(lat_subset),length(lon_subset));
                gsi_subset = NaN(n_years,length(lat_subset),length(lon_subset));
                vpd_subset = ones(length(lat_subset),1);

                % Call function, iterating over each lat/lon and year.
                for i=1:length(lon_subset)
                    for j=1:length(lat_subset)
                        for yr=1:n_years

                            lsf_subset(yr,j,i) = findlsf(m.data(:,yr,j,i));
                            gsi_subset(yr,j,i) = calcgsi(m.data(:,yr,j,i),...
                                                         day_subset(j,:),...
                                                         vpd_subset(j,:));
                        
                        end     % yr; n_years.
                    end         % j; lat_subset.
                end             % i; lon_subset.

                % Patch regional subsets together.
                lsf_CONUS(:,lat_subset,lon_subset) = lsf_subset;
                gsi_CONUS(:,lat_subset,lon_subset) = gsi_subset;

            end     % y; LAT_START.
        end         % x; LON_START.

        % Clear variables.
        clear model file_name m n_years lat x lon_subset y lat_subset day_subset 
        clear lsf_subset gsi_subset vpd_subset i j yr


	end 	        % e; Experiment loop.
end 		        % m; Model loop.










