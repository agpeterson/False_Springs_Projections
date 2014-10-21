%%=============================================================================
% NAME:   Process_GRIDMET.m
% AUTHOR: Alexander Peterson
% DATE:   20 Oct. 2014
% DESCR:  This script processes the GRIDMET data on Thunder to find last spring
%         freezes and GSI.
% IN:     GRIDMET
% OUT:    N/A
% CALLS:  findlsf.m; 
%==============================================================================

% Create path suffix and prefix strings to be concatenated.
PATH = '/storage/OBSDATA/METDATA/'
FILE = 'metdata_tmin_19792012_CONUS.mat'
WRITE_DIR = '/home/alex/';

% Break CONUS grid into regional subsets. The more powerful the machine, the
% larger the regional subsets can be.
LAT_START = [1 301];
LAT_END = [300 585];
LON_START = [1 301 601 901 1101];
LON_END = [300 600 900 1100 1386];

% Create constant for number of years, lat, lon, and models.
N_LAT = 585;
N_LON = 1386;
N_YRS = 34;
N_DAY = 366;

% Create year index.
YR_INDEX = 1:34;

%%=============================================================================
% Load or create necessary files.
%==============================================================================

% Create matfile pointers to files. 
lsf = matfile([WRITE_DIR,'gridMET_lsf.mat'],'Writable',true);
gsi = matfile([WRITE_DIR,'gridMET_gsi.mat'],'Writable',true);

% Preallocate space with NaNs.
lsf.lsf_CONUS = NaN(N_YRS,N_LAT,N_LON,'single');
gsi.gsi_CONUS = NaN(N_YRS,N_LAT,N_LON,'single');


%%=============================================================================
% Body of script to access and process MACA data. To do so, iterate over models
% and experiments to concatenate strings and access the .mat files, using
% parallel processing to subset CONUS and process the data.
%==============================================================================

% Open parallel processor pool.
if matlabpool('size') == 0

    matlabpool open local 12

end

% Create path string for each file and set pointer to matfile.
file_name = [PATH,FILE];
file = matfile(file_name);

% Create variable to hold number of years, switching on experiment.
n_yrs = 34;

% Create/load day_length variable (photoperiod) for use in GSI calculation. If
% file exists, load from file. If file does not exist, pull lat from one model
% and call calcdaylength function.

% Iterate over latitude to calculate day length for all days.
lat = file.lat;
lon = file.lon;

for i=1:N_LAT

    day_length(i,:) = calcdaylength(1:N_DAY,lat(i));

end
day_length = double(day_length');      % Transpose such that lat is outside.


%%=====================================================================
% Iterate over CONUS by breaking lat/lon into regional subset, then
% iterating over regional subset to call findlsf function.
%======================================================================

yr_index = [1:34];

for x=1:length(LON_START)
        
    % Break CONUS into regional lon subset.
    lon_subset = [LON_START(x):LON_END(x)];
        
    for y=1:length(LAT_START)
            
        % Write lon cell, lat cell, experiment, and model to output.
        [x y]
    
        % Break CONUS into regional lat, day, and vpd subsets.
        lat_subset = [LAT_START(y):LAT_END(y)];
        day_subset = day_length(:,lat_subset);
        vpd_subset = ones(N_DAY,length(lat_subset));    % Temporary.

        % Create temporary variable to store daily and yearly data for
        % each lat/lon subset.
       	t_var = double(file.data(:,yr_index,lat_subset,lon_subset));

        %%=============================================================
        % Process lon_subset using parallel function. To do so,
        % preallocate for LSF and GSI, then begin parfor loop iterating
        % over each lon_subset. Call findlsf and calcgsi on each
        % lon_subset, then concatenate together.
        %==============================================================
		
        % Create variables to store lengths of lat and lon.
        n_lat = length(lat_subset);
        n_lon = length(lon_subset);

        % Preallocate subset variables.
        lsf_subset = NaN(n_yrs,n_lat,n_lon);
        gsi_subset = NaN(n_yrs,n_lat,n_lon);
                
        tic
        % Parallel iteration over lon_subset.
        parfor i=1:n_lon
                    
            % Create temporary variable for each lon_subset, each
            % holding all latitudes for one longitude.
            t_var2 = t_var(:,:,:,i);

            % Call parallel function.
            [lsf_sub,gsi_sub] = findspringevents(t_var2,day_subset,vpd_subset);
                        
            % Concatenate subregional variables to subset.
            lsf_subset(:,:,i) = lsf_sub;
            gsi_subset(:,:,i) = gsi_sub;

        end     % i; lon_subset.
        toc

        % Write regional subsets together for all models.
        lsf.lsf_CONUS(yr_index,lat_subset,lon_subset) = single(lsf_subset);
        gsi.gsi_CONUS(yr_index,lat_subset,lon_subset) = single(gsi_subset);

    end     % y; LAT_START.
end         % x; LON_START.
    

matlabpool close;   % Close processor pool.


%% Extract 1979 to 2010 and save.
gridmet_gsi = gsi.gsi_CONUS(1:30,:,:);
gridmet_lsf = lsf.lsf_CONUS(1:30,:,:);
save('gridmet_19792009.mat','gridmet_gsi','gridmet_lsf')

