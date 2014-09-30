%%=============================================================================
% NAME:   Protoscript_v4.m
% AUTHOR: Alexander Peterson
% DATE:   21 Sept. 2014
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
MDL_TARGET = [1:20];
EXP_TARGET = [2:3];
VAR_TARGET = [2];

% Create path suffix and prefix strings to be concatenated.
PATH_PREFIX = '/storage/DOWNSCALED/CMIP5/MACAv2-METDATA/';
FILE_PREFIX = 'maca_v2_metdata_2var_10pat_CONUS_';
FILE_SUFFIX = {'_historical_tasmin.mat';...
               '_rcp45_tasmin.mat';...
               '_rcp85_tasmin.mat'};
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
N_MDL = 20;
N_YRS = 244;                   % Historical is 56 years, RCPs are both 94 years.
N_DAY = 365;

% Create index for historical and future years.
HST_INDEX = 1:56;
FUT_INDEX = [57:150] - 57 + 1;         % Sets FUT_INDEX as 1:94 for file access.


%%=============================================================================
% Load or create necessary files.
%==============================================================================

% Create/load day_length variable (photoperiod) for use in GSI calculation. If
% file exists, load from file. If file does not exist, pull lat from one model
% and call calcdaylength function.
if exist('day_length.mat') == 2

    load('day_length.mat','day_length')

else

    % Create file name string.
    model = char(MDL_NAME(MDL_TARGET(1)));
    file_name = [PATH_PREFIX,FILE_PREFIX,model,char(FILE_SUFFIX(1))];

    % Create file pointer and pull latitude.
    file = matfile(file_name);
    lat = file.lat;
    lon = file.lon;

    % Iterate over latitude to calculate day length for all days.
    for i=1:N_LAT
        day_length(i,:) = calcdaylength(1:N_DAY,lat(i));
    end
    day_length = double(day_length');      % Transpose such that lat is outside.

    % Save and clear variables.
    save('day_length.mat','day_length')
    clear model file_name file

end

% Create/load lsf and gsi matfiles. If lsf.mat and gsi.mat exist, create
% matfile pointers, else create new matlab files and fill with NaNs.
if exist('lsf.mat') == 2 && exist('gsi.mat') == 2

    lsf = matfile([WRITE_DIR,'lsf.mat'],'Writable',true);
    gsi = matfile([WRITE_DIR,'gsi.mat'],'Writable',true);

else

    % Preallocate space with NaNs.
    lsf.lsf_CONUS = NaN(N_YRS,N_LAT,N_LON,N_MDL,'single');
    gsi.gsi_CONUS = NaN(N_YRS,N_LAT,N_LON,N_MDL,'single');

    % Add model names to matfiles.
    lsf.MDL_NAME = MDL_NAME;
    gsi.MDL_NAME = MDL_NAME;

    % Add lat/lon to matfiles.
    lsf.lat = lat;
    lsf.lon = lon;
    gsi.lat = lat;
    gsi.lon = lon;

end


%%=============================================================================
% Body of script to access and process MACA data. To do so, iterate over models
% and experiments to concatenate strings and access the .mat files, using
% parallel processing to subset CONUS and process the data.
%==============================================================================

% Open parallel processor pool.
if matlabpool('size') == 0
    matlabpool open local 12
end

% Model and experiment iteration; run prototype run on CNRM-CM5.
m = 6   % for m=1:length(MDL_TARGET)
    
    % Subset model name using character array.
    model = char(MDL_NAME(MDL_TARGET(m)));
        
    for e = EXP_TARGET  % Possible break point?		
        
        % Create path string for each file and set pointer to matfile.
        file_name = [PATH_PREFIX,FILE_PREFIX,model,char(FILE_SUFFIX(e))];
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
            
                % Write lon cell, lat cell, experiment, and model to output.
                [x y e m]
    
                % Break CONUS into regional lat, day, and vpd subsets.
                lat_subset = [LAT_START(y):LAT_END(y)];
                day_subset = day_length(:,lat_subset);
                vpd_subset = ones(N_DAY,length(lat_subset));    % Temporary.

                % Create temporary variable to store daily and yearly data for
                % each lat/lon subset.
                if e == 1
                	t_var = double(file.data(:,HST_INDEX,...
                                               lat_subset,...
                                               lon_subset));
                else
                	t_var = double(file.data(:,FUT_INDEX,...
                                               lat_subset,...
                                               lon_subset));
                end

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
                    [lsf_sub,gsi_sub] = findfscomponents(t_var2,...
                                                         day_subset,...
                                                         vpd_subset);
                    
                    % Concatenate subregional variables to subset.
                    lsf_subset(:,:,i) = lsf_sub;
                    gsi_subset(:,:,i) = gsi_sub;

                end     % i; lon_subset.
                toc

                % Write regional subsets together for all models.
                lsf.lsf_CONUS(yr_index,lat_subset,lon_subset,m) = ...
                                                            single(lsf_subset);
                gsi.gsi_CONUS(yr_index,lat_subset,lon_subset,m) = ...
                                                            single(gsi_subset);

            end     % y; LAT_START.
        end         % x; LON_START.
     end            % e; Experiment loop.
 end 		        % m; Model loop.

matlabpool close;   % Close processor pool.


%%=============================================================================
% Find false springs using Spectrelight, not Thunder.
%==============================================================================

% Create constants, indices, though technically no need to separate based on
% RCP.
YRS_HST = 1:56;
YRS_R45 = 57:150;       % 57+94-1
YRS_R85 = 151:244;      % end-94+1
N_LAT = 585;
N_LON = 1386;
N_MDL = 20;
N_YRS = 244;            % Historical is 56 years, RCPs are both 94 years.
N_DAY = 365;

% Load GSI and LSF using matfile pointers.
file = matfile('gsi.mat');
gsi = file.gsi_CONUS(:,:,:,6);
clear file

file = matfile('lsf.mat');
lsf = file.lsf_CONUS(:,:,:,6);
clear file

% Find years with false springs across all lat/lon points.
fs = NaN(N_YRS,N_LAT,N_LON);
for i=1:N_LON
    for j=1:N_LAT
        for k=1:N_YRS
            if isnan(lsf(k,j,i)) || isnan(gsi(k,j,i))
                fs(k,j,i) = NaN;
            elseif lsf(k,j,i) >= (7+gsi(k,j,i))
                fs(k,j,i) = 1;
            else
                fs(k,j,i) = 0;
            end
        end
    end
end

% Separate based on RCP.
gsi_hst = gsi(YRS_HST,:,:);
gsi_r45 = gsi(YRS_R45,:,:);
gsi_r85 = gsi(YRS_R85,:,:);

lsf_hst = lsf(YRS_HST,:,:);
lsf_r45 = lsf(YRS_R45,:,:);
lsf_r85 = lsf(YRS_R85,:,:);

fs_hst = fs(YRS_HST,:,:);
fs_r45 = fs(YRS_R45,:,:);
fs_r85 = fs(YRS_R85,:,:);

save('gsi_v2.mat','gsi_hst','gsi_r45','gsi_r85')
save('lsf_v2.mat','lsf_hst','lsf_r45','lsf_r85')
save('fs_v2.mat','fs_hst','fs_r45','fs_r85')

% Calculate FSEI for 1950-2005 period.
fs_hst_sum = squeeze(sum(fs_hst,1));
lsf_hst_sum = squeeze(sum(~isnan(lsf_hst),1));
fsei_hst = (fs_hst_sum ./ lsf_hst_sum) * 100;
save('fs_v2.mat','fsei_hst','-append')

% Plot 1950-2005 GSI, LSF, and FSEI.
map_title = 'Mean Last Spring Freeze Date 1950-2005';
cb_type =
cb_color =
cb_units = 'Day of Year';
cb_flip = 'No Flip';
plotfscomponents(lsf_hst,1,181,20,colorbar_units,map_title,colorbar_flip)



% Map historical GSI.
gsi_hst = squeeze(nanmean(gsi,1));
map_title = 'Mean Plant Green-up Date 1950-2005';
colorbar_flip = 'No Flip';
plotfscomponents(gsi_hst,1,181,20,colorbar_units,map_title,colorbar_flip)

cb_type,cb_color,cb_units,cb_flip