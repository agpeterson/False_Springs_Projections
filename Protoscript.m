%%=============================================================================
% NAME:   Protoscript_v3.m
% AUTHOR: Alexander Peterson
% DATE:   16 Sept. 2014
% DESCR:  This script contains prototype code to read through one
%         MACAv2-METDATA model and calculate last spring freeze and GSI days.
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
%MDL_TARGET = [1:20];
MDL_TARGET = [6];       % CNRM-CM5 for prototype run.
EXP_TARGET = [1:2];
VAR_TARGET = [2];

% Create path suffix and prefix strings to be concatenated.
PATH_PREFIX = '/storage/DOWNSCALED/CMIP5/MACAv2-METDATA/';
FILE_PREFIX = 'maca_v2_metdata_2var_10pat_CONUS_';
FILE_SUFFIX = {'_historical_tasmin.mat';...
               '_rcp45_tasmin.mat';...
               '_rcp85_tasmin.mat'};

% Break CONUS grid into regional subsets. The more powerful the machine, the
% larger the regional subsets can be.
LAT_START = [1 201 401];
LAT_END = [200 400 585];
LON_START = [1 301 601 901 1101];
LON_END = [300 600 900 1100 1386];

% Create constant for number of years, lat, lon, and models.
N_LAT = 585;
N_LON = 1386;
N_MDL = 20;
N_YRS = 244;      % Historical is 56 years, RCPs are both 94 years.
N_DAY = 365;

% Load day_length variable.
load('day_length.mat','day_length')


%%=============================================================================
% Body of script to access and process MACA data. To do so, iterate over models
% and experiments to concatenate strings and access the .mat files, using
% parallel processing to subset CONUS and process the data.
%==============================================================================

% Open parallel processor pool.
if matlabpool('size') == 0
    matlabpool open local 12
end

% Create new matlab file to store CONUS variables.
write_dir = '/eddy/FALSE_SPRINGS/';
m_new = matfile([write_dir,'false_springs.mat'],'Writable',true);
m_new.lsf_CONUS = NaN(N_YRS,N_LAT,N_LON,N_MDL,'single');
m_new.gsi_CONUS = NaN(N_YRS,N_LAT,N_LON,N_MDL,'single');

% Model and experiment iteration. 
m = 1;  % for m=1:length(MDL_TARGET)

    % Subset model name using character array.
    model = char(MDL_NAME(MDL_TARGET(m)));
        
    e = 1; % for e=1:3  % REFACTOR - should not be hard-coded.
		
        % Create path string for each file.
        file_name = [PATH_PREFIX,FILE_PREFIX,model,char(FILE_SUFFIX(e))];
        		
        % Set pointer as matfile.
        file = matfile(file_name);

        % Create variable to hold number of years, switching on experiment.
        if e == 1
            yr_index = [1:56];
        elseif e == 2
            yr_index = [57:150];
        else e == 3
            yr_index = [151:244];
     	end;
        n_yrs = length(yr_index);

        %%=====================================================================
        % Iterate over CONUS by breaking lat/lon into regional subset, then
        % iterating over regional subset to call findlsf function.
        %======================================================================
        	
        for x=1:length(LON_START)
        
            % Break CONUS into regional lon subset.
            lon_subset = [LON_START(x):LON_END(x)];
        
            for y=1:length(LAT_START)
                
                % Break CONUS into regional lat, day, and vpd subsets.
                lat_subset = [LAT_START(y):LAT_END(y)];
                day_subset = day_length(:,lat_subset);
                vpd_subset = ones(N_DAY,length(lat_subset));

                % Call function, iterating over each lat/lon and year.
                t_var = double(file.data(:,:,lat_subset,lon_subset));

                %%=============================================================
                % Process lon_subset using parallel function. To do so,
                % preallocate for LSF and GSI, then begin parfor loop iterating
                % over each lon_subset. Call findlsf and calcgsi on each
                % lon_subset, then concatenate together.
                %==============================================================
		
                % Create variables to store lengths of lat and lon.
                n_lat=length(lat_subset);
                n_lon=length(lon_subset);

                % Preallocate subset variables.
                lsf_subset = NaN(n_yrs,n_lat,n_lon);
                gsi_subset = NaN(n_yrs,n_lat,n_lon);
                
                tic
                % Parallel iteration over lon_subset.
                parfor i=1:n_lon
                    
                    % Create temporary variable for each lon_subset.
                    t_var2 = t_var(:,:,:,i);

                    % Call parallel function.
                    [lsf_sub,gsi_sub] = fscomponents(t_var2,...
                                                     day_subset,...
                                                     vpd_subset);
                    
                    % Concatenate subregional variables to subset.
                    lsf_subset(:,:,i) = lsf_sub;
                    gsi_subset(:,:,i) = gsi_sub;

                end     % i; lon_subset.
                toc

                % Patch regional subsets together for all models.
                m_new.lsf_CONUS(yr_index,lat_subset,lon_subset,m) = ...
                                                            single(lsf_subset);
                m_new.gsi_CONUS(yr_index,lat_subset,lon_subset,m) = ...
                                                            single(gsi_subset);

            end     % y; LAT_START.
        end         % x; LON_START.
     end            % e; Experiment loop.
 end 		        % m; Model loop.

matlabpool close;   % Close processor pool.

