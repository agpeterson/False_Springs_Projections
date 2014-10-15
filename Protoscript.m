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
gsi = double(gsi);
clear file

file = matfile('lsf.mat');
lsf = file.lsf_CONUS(:,:,:,6);
lsf = double(lsf);
clear file


% Find years with false springs across all lat/lon points. --------------------
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

% Calculate FSEI for four periods: 1950-2005, 2010-2039, 2040-2069, 2070-2099.
YR_SERIES = horzcat([1950:2005],[2006:2099],[2006:2099]);
FUT_START = find(YR_SERIES == 2020);
FUT_END = find(YR_SERIES == 2049);
YRS_IND = {[1:56] [FUT_START(1):FUT_END(1)] [FUT_START(2):FUT_END(2)]};

% Preallocate, iterate over periods.
fs_sum = NaN(length(YRS_IND),N_LAT,N_LON);
fsei = NaN(length(YRS_IND),N_LAT,N_LON);
for i=1:length(YRS_IND)
    fs_sum(i,:,:) = sum(fs(YRS_IND{i},:,:),1);
    fsei(i,:,:) = (fs_sum(i,:,:) ./ sum(~isnan(lsf(YRS_IND{i},:,:)),1)) * 100;
end


% Find and plot normals and differences. --------------------------------------
% Calculate 1950-2005 mean.
lsf_hst_mean = squeeze(nanmean(lsf(YRS_HST,:,:),1));
gsi_hst_mean = squeeze(nanmean(gsi(YRS_HST,:,:),1));
fsei_hst = squeeze(fsei(1,:,:));

% Map historical values.
data = fsei_hst;
map_title = 'False Spring Exposure Index 1950-2005';
cb_type = 'seq';
cb_color = 'Blues';
cb_units = 'Percent of Years';
cb_flip = 'No Flip';
plotfscomponents(data,0,100,10,map_title,cb_type,cb_color,cb_units,cb_flip)


% Find difference between future and historical mean. -------------------------
lsf_diff = NaN(N_YRS,N_LAT,N_LON);
gsi_diff = NaN(N_YRS,N_LAT,N_LON);

YR_SERIES = horzcat([1950:2005],[2006:2099],[2006:2099]);
FUT_START = find(YR_SERIES == 2020);
FUT_END = find(YR_SERIES == 2049);
YRS_IND = {[FUT_START(1):FUT_END(1)] [FUT_START(2):FUT_END(2)]};

lsf_fut_diff = NaN(length(YRS_IND),N_LAT,N_LON);
gsi_fut_diff = NaN(length(YRS_IND),N_LAT,N_LON);

for i=1:length(YRS_IND)
    for j=1:N_YRS
        lsf_diff(j,:,:) = squeeze(lsf(j,:,:)) - lsf_hst_mean(:,:);
        gsi_diff(j,:,:) = squeeze(gsi(j,:,:)) - gsi_hst_mean(:,:);
    end
    lsf_fut_diff(i,:,:) = nanmean(lsf_diff(YRS_IND{i},:,:),1);
    gsi_fut_diff(i,:,:) = nanmean(gsi_diff(YRS_IND{i},:,:),1);
end

data = squeeze(lsf_fut_diff(1,:,:));
map_title = 'Difference in Mean Last Spring Freeze Date (2020-2049; RCP4.5)';
cb_type = 'div';
cb_color = 'RdBu';
cb_units = 'Departure from 1950-2005 Mean (Days)';
cb_flip = 'No Flip';
plotfscomponents(data,-30,30,6,map_title,cb_type,cb_color,cb_units,cb_flip)


% Relative change in FSEI between periods. New value minus old value, divided
% by old value.
fsei_r45 = squeeze(fsei(2,:,:));
fsei_r85 = squeeze(fsei(3,:,:));

fsei_r45_diff = fsei_r45 - fsei_hst;
fsei_r85_diff = fsei_r85 - fsei_hst;

% fsei_r45_change = ((fsei_r45 - fsei_hst) ./ fsei_hst)*100;
% fsei_r85_change = ((fsei_r85 - fsei_hst) ./ fsei_hst)*100;

% Map relative change.
data = fsei_r85_diff;
map_title = 'Difference in FSEI (2020-2049; RCP8.5)';
cb_type = 'div';
cb_color = 'RdBu';
cb_units = 'Departure from 1950-2005 FSEI (%)';
cb_flip = 'No Flip';
plotfscomponents(data,-50,50,10,map_title,cb_type,cb_color,cb_units,cb_flip)


% Create timeseries for Pullman. ----------------------------------------------
load conus_grid
pullman_lat = find(lat >= 46.7 & lat <= 46.75)
pullman_lon = find(lon >= -117.2 & lon <= -117.15)

pullman_gsi = gsi(:,pullman_lat,pullman_lon);
pullman_lsf = lsf(:,pullman_lat,pullman_lon);
pullman_fs = fs(:,pullman_lat,pullman_lon);

pullman_gsi_r45 = squeeze(pullman_gsi(1:150));
pullman_gsi_r85 = vertcat(pullman_gsi(1:56),pullman_gsi(151:end));

pullman_lsf_r45 = squeeze(pullman_lsf(1:150));
pullman_lsf_r85 = vertcat(pullman_lsf(1:56),pullman_lsf(151:end));


% Calculate moving 30-year sum.
SEMI_WINDOW = 15;
N_DATA = size(pullman_fs(1:150),1);  % size of data.
    
for i=1:150
    if i < SEMI_WINDOW + 1  % if i less than 16
        window = [1:i+SEMI_WINDOW]';    % window is 1:31
    elseif i > N_DATA - SEMI_WINDOW     % if i greater than 150 - 15
        window = [i-SEMI_WINDOW:N_DATA]';
    else
        window = [i-SEMI_WINDOW:i+SEMI_WINDOW]';
    end
    pullman_fsei_r45(i,:) = (sum(pullman_fs(window)) ./ length(window)) * 100;
end

% Plot.
SMOOTHING_WINDOW = 11;   % Number of years to smooth over; 5-1-5.
gsi_r45 = smooth(pullman_gsi_r45,SMOOTHING_WINDOW,'lowess');
gsi_r85 = smooth(pullman_gsi_r85,SMOOTHING_WINDOW,'lowess');
lsf_r45 = smooth(pullman_lsf_r45,SMOOTHING_WINDOW,'lowess');
lsf_r85 = smooth(pullman_lsf_r85,SMOOTHING_WINDOW,'lowess');
fsei_r45 = pullman_fsei_r45;
x = 1950:2099;


% Use plotyy to plot LSF and GSI on left y-axis, FSEI on right.
fig = figure('Position',[100 100 1000 618]);
set(gcf,'Color','W')
[ax,line1,line2] = ...
plotyy([x',x'],[gsi_r45,lsf_r45],x',fsei_r45,'plot','stairs')
axes(ax(1))
hold on
plot(x',gsi_r85,'linewidth',2,'linestyle','--','color',rgb('SteelBlue'));
plot(x',lsf_r85,'linewidth',2,'linestyle','--','color',rgb('IndianRed'));
set(ax(1),'ylim',[0 150],...
          'ytick',[0:30:150],...
          'ycolor','k',...
          'fontsize',12);
set(ax(2),'ylim',[0 100],...
          'ytick',[0:20:100],...
          'ycolor','k',...
          'fontsize',12);
set(line1,'linewidth',2);
set(line2,'linewidth',2,'color','k');
set(get(ax(1),'Ylabel'),'String','Day of Year','FontSize',14)
set(get(ax(2),'Ylabel'),'String','False Spring Exposure Index','FontSize',14)
h_legend = legend('RCP4.5 Green-up','RCP4.5 Last Spring Freeze',...
                  'RCP8.5 Green-up','RCP8.5 Last Spring Freeze','FSEI')
legend('boxoff')
set(h_legend,'FontSize',12);
title('Changes in LSF, GSI, and FS in Pullman, WA 1950-2099','FontSize',16)




% RCP8.5
fig = figure('Position',[100 100 1000 618]);
set(gcf,'Color','W')
[ax,line1,line2] = ...
plotyy([x',x'],[gsi_r85,lsf_r85],x',fsei_r45,'plot','stairs')
set(ax(1),'ylim',[0 150],...
          'ytick',[0:30:150],...
          'ycolor','k',...
          'fontsize',12);
set(ax(2),'ylim',[0 100],...
          'ytick',[0:20:100],...
          'ycolor','k',...
          'fontsize',12);
set(line1,'linewidth',2);
set(line2,'linewidth',2,'color','k');
set(get(ax(1),'Ylabel'),'String','Day of Year','FontSize',14)
set(get(ax(2),'Ylabel'),'String','False Spring Exposure Index','FontSize',14)
h_legend = legend('Green-up','Last Spring Freeze','FSEI')
legend('boxoff')
set(h_legend,'FontSize',12);
title('Bioclimatic Changes Under RCP8.5','FontSize',16)



% Check climatologies and tmin from MACA.
load pullman_tmin_C

% Convert to double.
tmin_hst = double(tmin_hst);
tmin_r45 = double(tmin_r45);
tmin_r85 = double(tmin_r85);

% Truncate data at 2099.
tmin_r45 = tmin_r45(:,1:end-1);
tmin_r85 = tmin_r85(:,1:end-1);

% Create 150 year time-series.
tmin_r45 = [tmin_hst,tmin_r45];
tmin_r85 = [tmin_hst,tmin_r85];

% Calculate annual tmins for JFMAM for years up to 2099.
hst_ann_mean = squeeze(nanmean(tmin_hst(1:151,:),1));
r45_ann_mean = squeeze(nanmean(tmin_r45(1:151,:),1));
r85_ann_mean = squeeze(nanmean(tmin_r85(1:151,:),1));

% Calculate climatological normal, take difference in annual tmins.
hst_normal = squeeze(nanmean(hst_ann_mean,2));
for i=1:150
    r45_ann_anom(i,:) = r45_ann_mean(:,i) - hst_normal;
    r85_ann_anom(i,:) = r85_ann_mean(:,i) - hst_normal;
end

% Plot timeseries for both absolutes and anomalies.
years = 1950:2099;
plot(years,r45_ann_anom)
hold all
plot(years,r85_ann_anom)
title('Changes in Mean JFMAM Tmin 1950-2099')
ylabel('Change in Degrees C')
grid on
legend('RCP4.5','RCP8.5','Location','Northwest')


% Daily anomalies - calculate normals for JFMAM days, then take difference.
daily_norm = squeeze(nanmean(tmin_hst(1:151,:),2));
for i=1:150
    r45_daily_anom(:,i) = tmin_r45(1:151,i) - daily_norm;
    r85_daily_anom(:,i) = tmin_r85(1:151,i) - daily_norm;
end

% Find coldest 10% of days each year and take mean.
for yr=1:150
    f1 = find(r45_daily_anom(:,yr) <= prctile(r45_daily_anom(:,yr),10,1));
    f2 = find(r85_daily_anom(:,yr) <= prctile(r85_daily_anom(:,yr),10,1));
    r45_extremes_anom(yr) = nanmean(r45_daily_anom(f1,yr));
    r85_extremes_anom(yr) = nanmean(r85_daily_anom(f2,yr));
end

plot(years,r45_extremes_anom)
hold all
plot(years,r85_extremes_anom)
title('JFMAM Coldest 10% of Days Anomalies 1950-2099')
xlabel('Year')
ylabel('Anomalies (Degrees C)')
grid on
legend('RCP4.5','RCP8.5','Location','Northwest')

scatter(r45_ann_anom,r45_extremes_anom)
title('RCP4.5 JFMAM Tmin Means vs Extremes 2006-2099')
xlabel('Mean Anomaly (Degrees C)')
ylabel('Coldest 10% Anomaly (Degrees C)')

scatter(r85_ann_anom,r85_extremes_anom)
title('RCP8.5 JFMAM Tmin Means vs Extremes 2006-2099')
xlabel('Mean Anomaly (Degrees C)')
ylabel('Coldest 10% Anomaly (Degrees C)')


% Difference between 2020-2049 and normal.
r45_extremes_normal = squeeze(nanmean(r45_extremes_anom(1:56)));
r85_extremes_normal = squeeze(nanmean(r85_extremes_anom(1:56)));

r45_extremes_diff = r45_extremes_anom - r45_extremes_normal;
r85_extremes_diff = r85_extremes_anom - r85_extremes_normal;

plot(years(71:100),r45_extremes_diff)
hold all
plot(years(71:100),r85_extremes_diff)
title('Anomalies in 10% Coldest Days')
ylabel('Difference from 1950-2005 Historical (Degrees C)')
legend('RCP45','RCP85','Location','Northwest')
grid on




% Look at USHCNv2 station data.
load('../False_Springs_Observed/Data_USHCN_18482013.mat')
load('../False_Springs_Observed/Data_StationInfo.mat')
years = 1848:2013;
find(years == 1950);
find(years == 2005);
tmin = tmindata(:,103:158,:);
vpd = ones(size(tmin,1),size(tmin,2),size(tmin,3));

% Daylight.
for i=1:1218
    daylit(i,:) = calcdaylength(1:366,lat(i));
end
daylit = repmat(daylit,[1 1 size(tmin,2)]);
daylit = permute(daylit,[2 3 1]);

% Find last spring freezes.
for i=1:1218
    for j=1:56
        [lsf(j,i)] = findlsf(tmin(:,j,i),-2.2);
    end
end

% GSI.
for i=1:1218
    gsiRaw(:,:,i) = calcgsi(tmin(:,:,i),daylit(:,:,i),vpd(:,:,i));
end
for i=1:1218
    gsi(:,i) = findgsi(gsiRaw(:,:,i));
end

% FS
for i=1:1218
    for j=1:56
        if isnan(lsf(j,i)) || isnan(gsi(j,i))
            fs(j,i) = NaN;
        elseif lsf(j,i) >= (7+gsi(j,i))
            fs(j,i) = 1;
        else
            fs(j,i) = 0;
        end
    end
end
for i=1:1218
    fs_sum(:,i) = sum(fs(:,i),1);
    fsei(:,i) = (fs_sum(:,i) ./ sum(~isnan(lsf(:,i)),1)) * 100;
end

s_fs = fs;
s_fs_sum = fs_sum;
s_fsei = fsei;
s_gsi = gsi;
s_lat = lat;
s_lon = lon;
s_lsf = lsf;
s_tmin = tmin;
s_vpd = vpd;

save('s_climatologies.mat','s_fs','s_fs_sum','s_fsei','s_gsi','s_lat',...
     's_lon','s_lsf','s_tmin','s_vpd')

% Find where stations and grid match; create variable with coordinate pairs,
% iterate over stations and find where lat/lon match within specified distance.
% If not found, input NaN.
load('s_climatologies.mat','s_lat','s_lon')
s_lat = s_lat';
s_lon = s_lon';

load('conus_grid.mat')
m_lat = lat;
m_lon = lon-360;

clear lat lon

stn_coords = [s_lon(:) s_lat(:)]
for i=1:1218
    p1 = find(m_lon >= stn_coords(i,1)-.03 & m_lon <= stn_coords(i,1)+.03);
    p2 = find(m_lat >= stn_coords(i,2)-.03 & m_lat <= stn_coords(i,2)+.03);
    if isempty(p1)
        p1 = find(m_lon == max(m_lon));
    elseif isempty(p2)
        p2 = find(m_lat == min(m_lat));
    end
    p3(i,:) = [p1(1) p2(1)];
end

lat_ind = round(p3(:,2));
lon_ind = round(p3(:,1));

m_lat2 = m_lat(lat_ind,:);
m_lon2 = m_lon(lon_ind,:);

figure();scatter(m_lat2,s_lat)
figure();scatter(m_lon2,s_lon)


% Pull GSI/LS/FS based on coordinate pairs.
file = matfile('gsi.mat');
gsi = file.gsi_CONUS(:,:,:,6);
gsi = double(gsi);
clear file
for i=1:1218
    m_gsi(i,:) = gsi(1:56,lat_ind(i),lon_ind(i));
end
clear gsi

file = matfile('lsf.mat');
lsf = file.lsf_CONUS(:,:,:,6);
lsf = double(lsf);
clear file
for i=1:1218
    m_lsf(i,:) = lsf(1:56,lat_ind(i),lon_ind(i));
end
clear lsf

for i=1:1218
    m_fs(i,:) = fs(1:56,lat_ind(i),lon_ind(i));
end
for i=1:1218
    fs_sum(:,i) = sum(m_fs(:,i),1);
    m_fsei(:,i) = (fs_sum(:,i) ./ sum(~isnan(m_lsf(:,i)),1)) * 100;
end



save('m_climatologies.mat','m_gsi','m_lsf')

% Create climatologies for comparison.
load('s_climatologies.mat','s_gsi','s_lsf')
load('m_climatologies.mat')
gsi_mean(1,:) = nanmean(s_gsi,1);
lsf_mean(1,:) = nanmean(s_lsf,1);
gsi_mean(2,:) = nanmean(m_gsi,2)';
lsf_mean(2,:) = nanmean(m_lsf,2)';

gsi_mean(:,155) = NaN;
lsf_mean(:,155) = NaN;

% Calculate statistics.
x1 = gsi_mean(1,:);
y1 = gsi_mean(2,:);
[gsi_r gsi_p] = corrcoef(x1,y1,'rows','pairwise')
p1 = polyfit(x1(~isnan(x1)),y1(~isnan(y1)),1)
yfit1 = polyval(p1,x1);
tbl = table(gsi_mean(1,:)',gsi_mean(2,:)','VariableNames',{'Obs','Mdl'});
mdl = fitlm(tbl)

x2 = lsf_mean(1,:);
y2 = lsf_mean(2,:);
[lsf_r lsf_p] = corrcoef(x2,y2,'rows','pairwise')
p2 = polyfit(x2(~isnan(x2)),y2(~isnan(y2)),1)
yfit2 = polyval(p2,x2);
tbl = table(lsf_mean(1,:)',lsf_mean(2,:)','VariableNames',{'Obs','Mdl'});
mdl = fitlm(tbl)

figure();
subplot(2,2,1:2)
plot(x1,y1,'+',x1,yfit1,'k','LineWidth',1.5)
title('Observed vs. Modeled GSI')
xlabel('Observed DoY'); ylabel('Modeled DoY')
legend('GSI','Regression','Location','Northwest')
grid on
text(140,65,{'RMSE: 3.64'; 'R2: 0.966'})

figure();
subplot(2,2,1:2)
plot(x2,y2,'+',x2,yfit2,'k','LineWidth',1.5)
title('Observed vs. Modeled LSF')
xlabel('Observed DoY'); ylabel('Modeled DoY')
legend('LSF','Regression','Location','Northwest')
grid on
text(140,30,{'RMSE: 15.8'; 'R2: 0.752'})

gsi_diff = gsi_mean(2,:) - gsi_mean(1,:);
lsf_diff = lsf_mean(2,:) - lsf_mean(1,:);

subplot(2,2,3)
plot(s_lat,gsi_diff,'+')
hold all
plot(25:50,zeros(1,26),'k','LineWidth',1.5)
title('Latitude vs. Difference between GSI')
xlabel('Latitude'); ylabel('Difference (Days)')
grid on

subplot(2,2,3)
plot(s_lat,lsf_diff,'+')
hold all
plot(25:50,zeros(1,26),'k','LineWidth',1.5)
title('Latitude vs. Difference between LSF')
xlabel('Latitude'); ylabel('Difference (Days)')
grid on

subplot(2,2,4)
plot(s_lon,gsi_diff,'+')
hold all
plot(-130:-60,zeros(1,71),'k','LineWidth',1.5)
title('Longitude vs. Difference between GSI')
xlabel('Longitude'); ylabel('Difference (Days)')
grid on

subplot(2,2,4)
plot(s_lon,lsf_diff,'+')
hold all
plot(-130:-60,zeros(1,71),'k','LineWidth',1.5)
title('Longitude vs. Difference between GSI')
xlabel('Longitude'); ylabel('Difference (Days)')
grid on

% Map relative change.
data = gsi_diff;
data_sig = ones(size(gsi_diff,1),size(gsi_diff,2));

plotstationdata(data,data_sig,lat,lon,30,-30);
title('Modeled GSI - Observed GSI')

data = lsf_diff;


x1 = s_fsei;
y1 = m_fsei;
[r p] = corrcoef(x1,y1,'rows','pairwise')
p1 = polyfit(x1(~isnan(x1)),y1(~isnan(y1)),1)
yfit1 = polyval(p1,x1);
tbl = table(x1,y1,'VariableNames',{'Obs','Mdl'});
mdl = fitlm(tbl)

validdata1 = ~isnan(x1)
validdata2 = ~isnan(y1)
validdataBoth = validdata1 & validdata2
keep1 = x1(validdataBoth) 
keep2 = y1(validdataBoth) 

p1 = polyfit(keep1, keep2, 1)
yfit1 = polyval(p1,x1);
tbl = table(x1',y1','VariableNames',{'Obs','Mdl'});
mdl = fitlm(tbl)


figure();
subplot(2,2,1:2)
plot(x1,y1,'+',x1,yfit1,'k','LineWidth',1.5)
title('Observed vs. Modeled FSEI')
xlabel('Observed %'); ylabel('Modeled %')
legend('FSEI','Regression','Location','Northwest')
grid on
text(80,10,{'RMSE: 12.5'; 'R2: 0.645'})

fsei_diff = m_fsei - s_fsei;

subplot(2,2,3)
plot(s_lat,fsei_diff,'+')
hold all
plot(25:50,zeros(1,26),'k','LineWidth',1.5)
title('Latitude vs. Difference between FSEI')
xlabel('Latitude'); ylabel('Difference (%)')
grid on

subplot(2,2,4)
plot(s_lon,fsei_diff,'+')
hold all
plot(-130:-60,zeros(1,71),'k','LineWidth',1.5)
title('Longitude vs. Difference between FSEI')
xlabel('Longitude'); ylabel('Difference (%)')
grid on


data = fsei_diff;
data_sig = ones(size(fsei_diff,1),size(fsei_diff,2));

plotstationdata(data,data_sig,s_lat,s_lon,60,-60);
title('Modeled FSEI - Observed FSEI')
















