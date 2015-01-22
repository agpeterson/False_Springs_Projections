%%=============================================================================
% Calculate VPD. Need specific humidity and pressure for calcDewPoint.
% Need tmax, tmin, tdew, and el for calcVPD.
% METDATA elevation located at /storage/OBSDATA/METDATA/dem/metdata_elev.nc
% tmax @ metdata_tmax_19792012_CONUS.mat
% tmin @ metdata_tmin_19792012_CONUS.mat
% sph @ metdata_sph_19792012_CONUS.mat
%==============================================================================
%%=============================================================================
% NAME:   Process_GRIDMET.m
% AUTHOR: Alexander Peterson
% DATE:   15 Nov. 2014
%
% DESCR:  This script processes the GridMET 1979-2012 data to find last spring
%         freezes and green-up.
% REF:    N/A
%
% IN:     metdata_tmin_19792012_CONUS.mat
%         metdata_tmax_19792012_CONUS.mat
%         metdata_sph_19792012_CONUS.mat
% OUT:    gridMET_vpd_gu.mat
% CALLS:  findSpringEvents.m (GridMET VPD specific)
%==============================================================================

% Break CONUS grid into regional subsets. The more powerful the machine, the
% larger the regional subsets can be.
LAT_START = [1 301];
LAT_END = [300 585];
LON_START = [1 601 1101];
LON_END = [600 1100 1386];

% Create constant for number of years, lat, lon, and models.
N_LAT = 585;
N_LON = 1386;
N_YRS = 34;
N_DAYS = 181;

% Create time indices.
YR_IND = 1:34;
DAY_IND = 1:181;


%%=============================================================================
% Access and create matfiles.
%==============================================================================

% Open GridMET matfiles.
sph = matfile('/storage/OBSDATA/METDATA/metdata_sph_19792012_CONUS.mat')
tmin = matfile('/storage/OBSDATA/METDATA/metdata_tmin_19792012_CONUS.mat')
tmax = matfile('/storage/OBSDATA/METDATA/metdata_tmax_19792012_CONUS.mat')

% Create files to store lsf and gsi.
gu = matfile('/home/alex/gridMET_vpd_gu.mat','Writable',true);

% Preallocate space with NaNs.
gu.gu_CONUS = NaN(N_YRS,N_LAT,N_LON,'single');

% Save lat and lon.
lat = tmin.lat;
lon = tmin.lon - 360;


%%=============================================================================
% Pressure.
%==============================================================================

% Open elevation data.
source = '/storage/OBSDATA/METDATA/dem/metdata_elev.nc'
ncid = netcdf.open(source,'NOWRITE')
ncdisp(source)
elv = ncread(source,'elevation');

% Calculate pressure.
pres = 1013.25 * (1-.0001 * elv);
pres = pres';


%%=============================================================================
% Photoperiod.
%==============================================================================

% Calculate photoperiod on tmin lat.
% for i=1:N_LAT
%     day_length(:,i) = calcDayLength(DAY_IND,lat(i));
% end
% save('day_length.mat','day_length')

% Open day_length data.
load('day_length.mat')

% Find where photoperiod is greater than or less than bounds and normalize; 
% all values below bounds set to 0 and all values above bounds set to 1.
HOUR_LOW = 10;
HOUR_HIGH = 11;
day_length = day_length(DAY_INDEX,:);
dayl = NaN(size(day_length)); 
f = find(day_length > HOUR_LOW & day_length < HOUR_HIGH);
dayl(f) = day_length(f) - HOUR_LOW;
dayl(day_length <= HOUR_LOW) = 0;
dayl(day_length >= HOUR_HIGH) = 1;


%%=============================================================================
% Body of script to access and process data.
%==============================================================================

% Open parallel processor pool.
if matlabpool('size') == 0

    matlabpool open local 12

end

% Iterate through regions.
for x=1:length(LON_START)
        
    % Break CONUS into regional lon subset.
    lon_subset = [LON_START(x):LON_END(x)];
    n_lon = length(lon_subset);

    for y=1:length(LAT_START)

        % Write lon and lat cell to output.
        [x y]
    
        % Break CONUS into regional lat subsets.
        lat_subset = [LAT_START(y):LAT_END(y)];
        n_lat = length(lat_subset);

        % Preallocate subset variables.
        gu_subset = NaN(N_YRS,n_lat,n_lon,'single');

        % Break pressure and day length into regional subset.
        pres_subset = pres(lat_subset,lon_subset);
        dayl_subset = dayl(:,lat_subset);

        % Create temporary variable to store daily and yearly data for
        % each lat/lon subset.
        t_sph = sph.data(DAY_IND,:,lat_subset,lon_subset);
        t_tmin = tmin.data(DAY_IND,:,lat_subset,lon_subset);
        t_tmax = tmax.data(DAY_IND,:,lat_subset,lon_subset);
        
        % Parallel iteration over lon_subset.
        tic
        parfor lon=1:n_lon

            % Call parallel function.
            [gu_sub] = findSpringEventsVPD(t_tmin(:,:,:,lon),...
                                           t_tmax(:,:,:,lon),...
                                           t_sph(:,:,:,lon),...
                                           pres_subset(:,lon),...
                                           dayl_subset(:,:),...
                                           N_DAYS,N_YRS,n_lat,1);

            % Concatenate subregional variables to subset.
            gu_subset(:,:,lon) = gu_sub;

        end     % parfor - lon; n_lon.
        toc

        % Write regional subsets to file.
        gu.gu_CONUS(YR_IND,lat_subset,lon_subset) = gu_subset;

    end             % y; LAT_START.
end                 % x; LON_START.

matlabpool close;   % Close processor pool.