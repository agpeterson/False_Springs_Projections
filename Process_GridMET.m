%%=============================================================================
% NAME:   Process_GRIDMET.m
% AUTHOR: Alexander Peterson
% DATE:   15 Nov. 2014
% DESCR:  This script processes the GRIDMET data on Thunder to find last spring
%         freezes and green-up.
% IN:     GRIDMET
% OUT:    N/A
% CALLS:  
%==============================================================================

% Create path suffix and prefix strings to be concatenated.
PATH = '/storage/OBSDATA/METDATA/'
FILE = 'metdata_tmin_19792012_CONUS.mat'
WRITE_DIR = '/home/alex/';

% Break CONUS grid into regional subsets. The more powerful the machine, the
% larger the regional subsets can be.
LAT_START = [1];
LAT_END = [585];
LON_START = [1 601 1101];
LON_END = [600 900 1386];

% Create constant for number of years, lat, lon, and models.
N_LAT = 585;
N_LON = 1386;
N_YRS = 34;
N_DAYS = 181;

% Create indices.
YR_INDEX = 1:34;
DAY_INDEX = 1:181;


%%=============================================================================
% Create lsf and gsi files.
%==============================================================================

% Create matfile pointers to files. 
lsf = matfile([WRITE_DIR,'gridMET_vpd_lsf.mat'],'Writable',true);
gsi = matfile([WRITE_DIR,'gridMET_vpd_gsi.mat'],'Writable',true);

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

% Iterate over latitude to calculate day length for all days.
lat = file.lat;
lon = file.lon;

for i=1:N_LAT

    day_length(i,:) = calcdaylength(1:N_DAYS,lat(i));

end
day_length = double(day_length');      % Transpose such that lat is outside.

% Photoperiod. Find where photoperiod is greater than or less than bounds and
% normalize; all values below bounds set to 0 and all values above bounds set
% to 1.
HOUR_LOW = 10;
HOUR_HIGH = 11;
day_length = day_length(DAY_INDEX,:);
dayl = NaN(size(day_length)); 
f = find(day_length > HOUR_LOW & day_length < HOUR_HIGH);
dayl(f) = day_length(f) - HOUR_LOW;
dayl(day_length <= HOUR_LOW) = 0;
dayl(day_length >= HOUR_HIGH) = 1;

% VPD. Assume held constant across time and space.
VPD_LOW = 0.9;
VPD_HIGH = 4.1;
vpd = single(1.0 - (1.0 - VPD_LOW) / (VPD_HIGH-VPD_LOW));

% Preallocate GSI_MAX.
GSI_MAX = NaN(N_LAT,N_LON,'single');

SCN = 1;


%%=====================================================================
% Iterate over CONUS by breaking lat/lon into regional subset, then
% iterating over regional subset to call findlsf function.
%======================================================================

for x=1:length(LON_START)
        
    % Break CONUS into regional lon subset.
    lon_subset = [LON_START(x):LON_END(x)];
    n_lon = length(lon_subset);
        
    for y=1:length(LAT_START)
            
        % Write lon cell, lat cell, experiment, and model to output.
        [x y]
    
        % Break CONUS into regional lat, day, and vpd subsets.
        lat_subset = [LAT_START(y):LAT_END(y)];
        n_lat = length(lat_subset);

        % Create temporary variable to store daily and yearly data for
        % each lat/lon subset.
       	t_var = file.data(DAY_INDEX,YR_INDEX,lat_subset,lon_subset);


        %%=============================================================
        % Process lon_subset using parallel function. To do so,
        % preallocate for LSF and GSI, then begin parfor loop iterating
        % over each lon_subset. Call findlsf and calcgsi on each
        % lon_subset, then concatenate together.
        %==============================================================

        % Preallocate subset variables.
        lsf_subset = NaN(N_YRS,n_lat,n_lon,'single');
        gsi_subset = NaN(N_YRS,n_lat,n_lon,'single');
        gsi_max_subset = NaN(n_lat,n_lon,'single');
                
        tic
        % Parallel iteration over lon_subset.
        parfor lon=1:n_lon
                    
            % Call parallel function.
            [lsf_sub,gsi_sub,gsi_max_sub] = findSpringEvents(t_var(:,:,:,lon),...
                GSI_MAX(lat_subset,lon_subset(lon)),dayl,vpd,N_DAYS,N_YRS,...,
                n_lat,YR_INDEX,SCN);
                        
            % Concatenate subregional variables to subset.
            lsf_subset(:,:,lon) = lsf_sub;
            gsi_subset(:,:,lon) = gsi_sub;
            gsi_max_subset(:,lon) = gsi_max_sub;

        end     % lon; lon_subset.
        toc

        % Write regional gsi_max subsets.
        GSI_MAX(lat_subset,lon_subset) = gsi_max_subset;

        % Write regional subsets together.
        lsf.lsf_CONUS(YR_INDEX,lat_subset,lon_subset) = lsf_subset;
        gsi.gsi_CONUS(YR_INDEX,lat_subset,lon_subset) = gsi_subset;

    end     % y; LAT_START.
end         % x; LON_START.
    

matlabpool close;   % Close processor pool.
