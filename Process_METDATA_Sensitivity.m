%%=============================================================================
% NAME:   Process_METDATA_Sensitivity.m
% AUTHOR: Alexander Peterson
% DATE:   15 Nov. 2014
%
% DESCR:  This script processes the GridMET 1979-2012 data to find last spring
%         freezes and green-up.
% REF:    N/A
%
% IN:     metdata_tmin_19792012_CONUS.mat
% OUT:    gridMET_lsf.mat; gridMET_gsi.mat
% CALLS:  findSpringEvents.m (GridMET specific)
%==============================================================================

% Break CONUS grid into regional subsets. The more powerful the machine, the
% larger the regional subsets can be.
LAT_START = [1];
LAT_END = [585];
LON_START = [1 701];
LON_END = [700 1386];

% Create constant for number of years, lat, lon, and models.
N_LAT = 585;
N_LON = 1386;
N_YRS = 30;
N_DAYS = 181;
N_MDL = 20;

% Create time indices.
YR_IND = 2:31;
DAY_IND = 1:181;


%%=============================================================================
% Access and create matfiles.
%==============================================================================

% Open GridMET matfiles.
disp('Accessing METDATA file...')
tmin_file = matfile('/storage/OBSDATA/METDATA/metdata_tmin_19792012_CONUS.mat')

% Create files to store lsf and gsi.
disp('Creating and preallocating GU and LSF files...')
lsf_file = matfile('/home/alex/LSF_Sensitivity.mat','Writable',true);
gu_file = matfile('/home/alex/GU_Sensitivity.mat','Writable',true);

% Preallocate space with NaNs.
lsf_file.data = NaN(N_YRS,N_LAT,N_LON,N_MDL,'single');
gu_file.data = NaN(N_YRS,N_LAT,N_LON,N_MDL,'single');

% Load tmin deltas for sensitivity experiment.
disp('Loading MACA tmin delta...')
load('Tmin_Delta.mat');


%%=============================================================================
% GSI constants.
%==============================================================================

disp('Creating day length and VPD variables...')

% Calculate photoperiod on tmin lat.
% for i=1:N_LAT
%     day_length(:,i) = calcDayLength(DAY_IND,lat(i));
% end
% save('day_length.mat','day_length')

% Open day_length data.
load('Day_Length_METDATA.mat')
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


%%=============================================================================
% Body of script to access and process data.
%==============================================================================

% Open parallel processor pool.
disp('Opening parallel processor pool...')
if matlabpool('size') == 0
    matlabpool open local 12
end

% Iterate through regions.
disp('Entering main loop...')
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
        lsf_subset = NaN(N_YRS,n_lat,n_lon,'single');
        gu_subset = NaN(N_YRS,n_lat,n_lon,'single');

        % Break day length into regional subset.
        dayl_subset = dayl(:,lat_subset);

        % Create temporary variable to store daily and yearly data for
        % each lat/lon subset. Replace -9999 (missing values) with NaN.
        disp('Loading temperature data...')
        t_tmin = tmin_file.data(DAY_IND,YR_IND,lat_subset,lon_subset);
        t_tmin(t_tmin == -9999) = NaN;

        % Model loop.
        disp('Entering model loop...')
        for mdl=1:N_MDL

            disp({'Model: ',mdl})

            % Extract tmin deltas for lat/lon subsets.
            t_delta = squeeze(tmin_delta(lat_subset,lon_subset,mdl));

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
                lsf_subset(:,:,lon) = lsf_sub;
                gu_subset(:,:,lon) = gu_sub;

            end     % parfor - lon; n_lon.
            toc

            % Write regional subsets to file.
            disp('Writing output to file...')
            lsf_file.data(1:30,lat_subset,lon_subset,mdl) = lsf_subset;
            gu_file.data(1:30,lat_subset,lon_subset,mdl) = gu_subset;

        end     % mdl; 1:N_MDL

    end         % y; LAT_START.
end             % x; LON_START.

% Close parallel pool and program.
matlabpool close;
exit