%%=====================================================================
% Calculate VPD. Need specific humidity and pressure for calcDewPoint.
% Need tmax, tmin, tdew, and el for calcVPD.
% METDATA elevation located at /storage/OBSDATA/METDATA/dem/metdata_elev.nc
% tmax @ metdata_tmax_19792012_CONUS.mat
% tmin @ metdata_tmin_19792012_CONUS.mat
% sph @ metdata_sph_19792012_CONUS.mat
%======================================================================

% Open GridMET matfiles.
sph = matfile('/storage/OBSDATA/METDATA/metdata_sph_19792012_CONUS.mat')
tmin = matfile('/storage/OBSDATA/METDATA/metdata_tmin_19792012_CONUS.mat')
tmax = matfile('/storage/OBSDATA/METDATA/metdata_tmax_19792012_CONUS.mat')

% Break CONUS grid into regional subsets. The more powerful the machine, the
% larger the regional subsets can be.
LAT_START = [1];
LAT_END = [585];
LON_START = [1 301 601 901 1101];
LON_END = [300 600 900 1100 1386];

% Create constant for number of years, lat, lon, and models.
N_LAT = 585;
N_LON = 1386;
N_YRS = 34;
N_DAYS = 181;

% Create indices.
YR_INDEX = 1:34;
DAY_INDEX = 1:181;

% Calculate tdew.

% Preallocate tdew.
tdew = matfile('/home/alex/gridMET_tdew.mat','Writable',true);

% Preallocate space with NaNs.
tdew = NaN(N_DAYS,N_YRS,N_LAT,N_LON,'single');
vpd = NaN(N_DAYS,N_YRS,N_LAT,N_LON,'single');

% Open elevation data.
source = '/storage/OBSDATA/METDATA/dem/metdata_elev.nc'
ncid = netcdf.open(source,'NOWRITE')
ncdisp(source)
elv = ncread(source,'elevation');

% Calculate pressure.
pres = 1013.25 * (1-.0001 * elv);
pres = pres';

% Iterate through regions.
for x=1:length(LON_START)
        
    % Break CONUS into regional lon subset.
    lon_subset = [LON_START(x):LON_END(x)];
    n_lon = length(lon_subset);

    % Create temporary variable to store daily and yearly data for
    % each lat/lon subset.
    t_sph = double(sph.data(DAY_INDEX,:,:,lon_subset));
    t_pres = pres(:,lon_subset);

    % Preallocate subset variables.
    tdew_subset = NaN(N_DAYS,N_YRS,N_LAT,n_lon,'single');

    tic
    % Parallel iteration over lon_subset.
    parfor lon=1:n_lon
        for lat=1:N_LAT
                    
            % Call parallel function for dew point, VPD, and spring events.
            [tdew_sub] = calcDewPoint(t_sph(:,:,lat,lon),t_pres(lat,lon));
                            
            % Concatenate subregional variables to subset.
            tdew_subset(:,:,lat,lon) = single(tdew_sub);

        end     % lat; lat_subset.
    end         % lon; lon_subset.
    toc

    % Write regional gsi_max subsets.
    tdew.tdew_CONUS(DAY_INDEX,:,:,lon_subset) = tdew_subset;

end
