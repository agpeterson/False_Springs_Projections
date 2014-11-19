%%=============================================================================
% NAME:   Process_MACAv2METDATA.m
% AUTHOR: Alexander Peterson
% DATE:   21 Oct. 2014
% DESCR:  This script reads through all downscaled MACAv2-METDATA models and 
%         calculates last spring freeze and green-up days.
% IN:     MACAv2-METDATA
% OUT:    gsi.mat; lsf.mat
% CALLS:  findlsf.m; findspringevents.m; findgreenup.m; calcgsi.m
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
scn_NAME = {'rcp45';...
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
MDL_TARGET = [1];
scn_TARGET = [2:3];
VAR_TARGET = [2];

% Create path suffix and prefix strings to be concatenated.
PATH_PREFIX = '/storage/DOWNSCALED/CMIP5/MACAv2-METDATA/';
FILE_PREFIX = 'maca_v2_metdata_2var_10pat_CONUS_';
FILE_SUFFIX = {'_historical_tasmin.mat';...
               '_rcp45_tasmin.mat';...
               '_rcp85_tasmin.mat'};
WRITE_DIR = '/home/katherine/ALEX/DATA/';

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
N_YRS = 244;                   % Historical is 56 years, RCPs are both 94 years.
N_DAYS = 181 

% Create index for historical and future years.
HST_INDEX = 1:56;
FUT_INDEX = [57:150] - 57 + 1;         % Sets FUT_INDEX as 1:94 for file access.
DAY_INDEX = 1:181;


%%=============================================================================
% Photoperiod.
%==============================================================================

% Load day_length variable (photoperiod) for use in GSI calculation. 
load('day_length.mat','day_length')


% Find where photoperiod is greater than or less than bounds and
% normalize; all values below bounds set to 0 and all values above bounds set
% to 1.
HOUR_LOW = 10;
HOUR_HIGH = 11;
day_length=day_length(DAYS,:);
dayl=NaN(size(day_length));%preallocate
f = find(day_length > HOUR_LOW &day_length < HOUR_HIGH);
dayl(f) = (day_length(f)-HOUR_LOW);
dayl(day_length <= HOUR_LOW) = 0;
dayl(day_length >= HOUR_HIGH) = 1;


%%=============================================================================
% VPD.
%==============================================================================
% Assume VPD is one for all lat/lon.
VPD_LOW = 0.9;
VPD_HIGH = 4.1;
vpd= single(1.0 -(1.0-VPD_LOW) / (VPD_HIGH-VPD_LOW));


%%=============================================================================
% Parallel processing.
%==============================================================================

% Set computer - either thunder or graupelX
COMPUTER = 'graupel1';

% Open parallel processor pool.
if strcmp(COMPUTER,'thunder')
	if matlabpool('size') == 0
        matlabpool open local 12
	end

elseif strcmp(COMPUTER,'graupel1')
	c = parcluster('local')
    c.NumWorkers = 12
    matlabpool open local 12

elseif strcmp(COMPUTER,'graupel2')
	c = parcluster('local2')
    c.NumWorkers = 12
    matlabpool open local2 12

elseif strcmp(COMPUTER,'graupel3')
	c = parcluster('local3')
    c.NumWorkers = 12
    matlabpool open local3 12

end


%%=============================================================================
% Body of script to access and process MACA data. To do so, iterate over models
% and experiments to concatenate strings and access the .mat files, using
% parallel processing to subset CONUS and process the data.
%==============================================================================

% Model iteration.
for m=MDL_TARGET

    % Create filename character strings.
	lsf_filename = [WRITE_DIR,'lsf','_',char(MDL_NAME{m}),'.mat'];
	gsi_filename = [WRITE_DIR,'gsi','_',char(MDL_NAME{m}),'.mat'];
	
    % Create/load lsf and gsi matfiles. If lsf.mat and gsi.mat exist, create
    % matfile pointers, else create new matlab files and fill with NaNs.
    if ~exist(lsf_filename) || ~exist(gsi_filename)
		lsf = matfile(lsf_filename,'Writable',true);
		gsi = matfile(gsi_filename,'Writable',true);

		% Preallocate space with NaNs.
		lsf.lsf_CONUS = NaN(N_YRS,N_LAT,N_LON,'single');
		gsi.gsi_CONUS = NaN(N_YRS,N_LAT,N_LON,'single');

		% Add lat/lon to matfiles.
		lsf.lat = lat;
		lsf.lon = lon - 360;
		gsi.lat = lat;
		gsi.lon = lon - 360;
	
    else
  		lsf = matfile(lsf_filename,'Writable',true);
        gsi = matfile(gsi_filename,'Writable',true);
	
    end

    % Preallocate GSI_MAXCLIM.
    GSI_MAXCLIM = NaN(N_LAT,N_LON,'single');
    
    % Subset model name using character array.
    model = char(MDL_NAME(m));
        
    % Switch on experiment.
    for scn=[1 length(FILE_SUFFIX)];
        
        % Create path string for each file and set pointer to matfile.
        file_name = [PATH_PREFIX,FILE_PREFIX,model,char(FILE_SUFFIX(scn))];
        file = matfile(file_name);

        % Create variable to hold number of years, switching on experiment.
        if scn == 1
            yr_index = [1:56];
		
        elseif scn == 2
            yr_index = [57:150];

        else scn == 3
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
            n_lon = length(lon_subset);
        
            for y=1:length(LAT_START)
            
                % Write lat/lon cell, experiment, and model to output.
                [m scn y x]
    
                % Break CONUS into regional lat, day, and vpd subsets.
                lat_subset = [LAT_START(y):LAT_END(y)];
                n_lat = length(lat_subset);

                % Create temporary variable to store daily and yearly data for
                % each lat/lon subset.
                if scn == 1
                	t_var = file.data(DAYS,HST_INDEX,lat_subset,lon_subset);
                
                else
                	t_var = file.data(DAYS,FUT_INDEX,lat_subset,lon_subset);
                end
                
                %%=============================================================
                % Process lon_subset using parallel function. To do so,
                % preallocate for LSF and GSI, then begin parfor loop iterating
                % over each lon_subset. Call findlsf and calcgsi on each
                % lon_subset, then concatenate together.
                %==============================================================
                
                % Preallocate subset variables.
                lsf_subset = NaN(n_yrs,n_lat,n_lon,'single');
                gsi_subset = NaN(n_yrs,n_lat,n_lon,'single');
                gsi_subset_maxclim = NaN(n_lat,n_lon,'single');
                
                tic 
                % Parallel iteration over lon_subset.
                parfor lon=1:n_lon

                    % Call parallel function.
                    [lsf_sub,...
                    gsi_sub,...
                    gsi_maxclim] = findSpringEvents(t_var(:,:,:,lon),...
                                    GSI_MAXCLIM(lat_subset,lon_subset(lon)),...
                                    dayl,vpd,n_days,n_yrs,n_lat,...
                                    HST_INDEX,scn);
                    
                    % Concatenate subregional variables to subset.
                    lsf_subset(:,:,lon) = lsf_sub;
                    gsi_subset(:,:,lon) = gsi_sub;
                    gsi_subset_maxclim(:,lon) = gsi_maxclim;

                end     % i; lon_subset.
                toc

                if scn==1

                    GSI_MAXCLIM(lat_subset,lon_subset) = gsi_subset_maxclim;

                end

                % Write regional subsets together for all models.
                lsf.lsf_CONUS(yr_index,lat_subset,lon_subset) = lsf_subset;
                gsi.gsi_CONUS(yr_index,lat_subset,lon_subset) = gsi_subset;

            end     % y; LAT_START.
        end         % x; LON_START.
    end             % e; Experiment loop.
end 		        % m; Model loop.

matlabpool close;   % Close processor pool.

