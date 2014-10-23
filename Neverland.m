% Neverland script



% Find differences between historical and future for RCP4.5.
% Load GSI and LSF using matfile pointers.
file = matfile('gsi.mat');
gsi = file.gsi_CONUS(57:150,:,:,6);
gsi = double(gsi);
clear file

file = matfile('lsf.mat');
lsf = file.lsf_CONUS(57:150,:,:,6);
lsf = double(lsf);
clear file

% Find difference between future years and 1950-2005 normal.
gsi_diff = NaN(94,N_LAT,N_LON);
lsf_diff = NaN(94,N_LAT,N_LON);

for i=1:94
    gsi_diff(i,:,:) = squeeze(gsi(i,:,:)) - gsi_hst;
    lsf_diff(i,:,:) = squeeze(lsf(i,:,:)) - lsf_hst;
end

% Average the difference across 3 futures: 2010-2039, 2040-2069, 2070-2099.
fut1 = [5:34];
fut2 = [35:64];
fut3 = [64:94];
futures = {fut1 fut2 fut3};

gsi_fut = NaN(3,N_LAT,N_LON);
lsf_fut = NaN(3,N_LAT,N_LON);
for i=1:3
    gsi_fut(i,:,:) = nanmean(gsi_diff(futures{i},:,:),1);
    lsf_fut(i,:,:) = nanmean(lsf_diff(futures{i},:,:),1);
end

% Plot differences.
lsf_fut3 = squeeze(lsf_fut(3,:,:));
map_title = 'Change in Last Spring Freeze (2070-2099; RCP4.5)';
cb_units = 'Change in Days';
cb_name = 'Blue';
cb_flip = 'Flip';
plotfscomponents(lsf_fut3,-80,0,10,map_title,cb_units,cb_name,cb_flip)

gsi_fut3 = squeeze(gsi_fut(3,:,:));
map_title = 'Change in Plant Green-up (2070-2099; RCP4.5)';
cb_units = 'Change in Days';
cb_name = 'Blue';
cb_flip = 'Flip';
plotfscomponents(gsi_fut3,-80,0,10,map_title,cb_units,cb_name,cb_flip)

% Calculate FSEI for each future.
% First find false springs.
for a=1:N_LON
    for b=1:N_LAT
        for c=1:94
            if isnan(lsf(c,b,a)) || isnan(gsi(c,b,a))
                fs(c,b,a) = NaN;
            elseif lsf(c,b,a) >= (7+gsi(c,b,a))
                fs(c,b,a) = 1;
            else
                fs(c,b,a) = 0;
            end
        end
    end
end

fs_sum = squeeze(sum(fs(futures{1},:,:),1));
fsei = (fs_sum./ 30)*100;

% Find difference in FSEI.
fsei_diff = fsei - fsei_hst;

map_title = 'Change in FSEI (2010-2039; RCP4.5)';
cb_units = 'Difference (Percent)';
cb_name = 'RedBlue';
cb_flip = 'Flip';
plotfscomponents(fsei_diff,-100,100,10,map_title,cb_units,cb_name,cb_flip)


gsi_r45 = double(gsi_r45);
gsi_r85 = double(gsi_r85);

lsf_r45 = double(lsf_r45);
lsf_r85 = double(lsf_r85);

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




fs_hst_sum = squeeze(sum(fs_hst,1));
lsf_hst_sum = squeeze(sum(~isnan(lsf_hst),1));
fsei_hst = (fs_hst_sum ./ lsf_hst_sum) * 100;
save('fs_v2.mat','fsei_hst','-append')



save('gsi.mat','gsi_diff','gsi_hst_mean','-append')
save('lsf.mat','lsf_diff','lsf_hst_mean','-append')
save('fs.mat','fs','fs_sum','fsei')

% Plot 1950-2005 GSI, LSF, and FSEI.
% LSF
data = squeeze(nanmean(lsf_hst,1));
map_title = 'Mean Last Spring Freeze Date 1950-2005';
cb_type = 'div';
cb_color = 'RdBu';
cb_units = 'Day of Year';
cb_flip = 'No Flip';
plotfscomponents(data,1,181,20,map_title,cb_type,cb_color,cb_units,cb_flip)

% GSI
data = squeeze(nanmean(gsi_hst,1));
map_title = 'Mean Plant Green-up Date 1950-2005';
cb_type = 'div';
cb_color = 'RdBu';
cb_units = 'Day of Year';
cb_flip = 'No Flip';
plotfscomponents(data,1,181,20,map_title,cb_type,cb_color,cb_units,cb_flip)

% FSEI
data = squeeze(fsei(1,:,:));
map_title = 'False Spring Exposure Index 1950-2005';
cb_type = 'div';
cb_color = 'RdBu';
cb_units = 'Percent of Years';
cb_flip = 'Flip';
plotfscomponents(data,0,100,10,map_title,cb_type,cb_color,cb_units,cb_flip)

% Calculate mean differences for each period.
YRS_IND = {[61:90] [91:120] [121:150] [155:184] [185:214] [215:244]};

gsi_diff_mean = NaN(7,N_LAT,N_LON);
lsf_diff_mean = NaN(7,N_LAT,N_LON);
for c=1:length(YRS_IND)
    gsi_diff_mean(c,:,:) = nanmean(gsi_diff(YRS_IND{c},:,:),1);
    lsf_diff_mean(c,:,:) = nanmean(lsf_diff(YRS_IND{c},:,:),1);
end

% Plot differences as boxplots? Calculate trends for maps?
% LSF 2010-2039.
data = squeeze(lsf_diff_mean(3,:,:));
map_title = 'Mean Last Spring Freeze Difference (2040-2070; RCP4.5)';
cb_type = 'div';
cb_color = 'RdBu';
cb_units = 'Difference';
cb_flip = 'No Flip';
plotfscomponents(data,-90,0,10,map_title,cb_type,cb_color,cb_units,cb_flip)

% GSI 2010-2039.
data = squeeze(gsi_diff_mean(3,:,:));
map_title = 'Mean Green-up Difference (2010-2039; RCP4.5)';
cb_type = 'div';
cb_color = 'RdBu';
cb_units = 'Difference';
cb_flip = 'No Flip';
plotfscomponents(data,-90,0,10,map_title,cb_type,cb_color,cb_units,cb_flip)


% Linear trend for plots; boxplots for periods.
gsi_r85_trend = NaN(N_LAT,N_LON);
gsi_r85_sig = NaN(N_LAT,N_LAT);
lsf_r85_trend = NaN(N_LAT,N_LON);
lsf_r85_sig = NaN(N_LAT,N_LON);
for i=1:N_LON
    for j=1:N_LAT
        [gsi_r85_trend(j,i),gsi_r85_sig(j,i)] = trenddata(gsi(155:244,j,i));
        [lsf_r85_trend(j,i),lsf_r85_sig(j,i)] = trenddata(lsf(155:244,j,i));
    end
end


data = gsi_r85_trend;
map_title = 'Plant Green-up Trend (2010-2099; RCP8.5)';
cb_type = 'seq';
cb_color = 'Reds';
cb_units = 'Trend (days over 90 years)';
cb_flip = 'Flip';
plotfscomponents(data*90,-90,0,10,map_title,cb_type,cb_color,cb_units,cb_flip)


data = lsf_r85_trend;
map_title = 'Last Spring Freeze Trend (2010-2099; RCP8.5)';
cb_type = 'seq';
cb_color = 'Reds';
cb_units = 'Trend (days over 90 years)';
cb_flip = 'Flip';
plotfscomponents(data*90,-90,0,10,map_title,cb_type,cb_color,cb_units,cb_flip)


% Relative difference. Take difference between all years, find climatological
% normal, then take future period differences.
rd = NaN(N_YRS,N_LAT,N_LON);
for i=1:N_LON
    for j=1:N_LAT
        rd(:,j,i) = lsf(:,j,i) - gsi(:,j,i);
    end
end

% Climatological normal.
rd_hst_mean = NaN(N_LAT,N_LON);
rd_diff = NaN(N_YRS,N_LAT,N_LON);
rd_hst_mean = squeeze(nanmean(rd(1:56,:,:),1));
for i=1:N_YRS
    rd_diff(i,:,:) = squeeze(rd(i,:,:)) - rd_hst_mean(:,:);
end

% Calculate mean differences for each period.
YRS_IND = {[61:90] [91:120] [121:150] [155:184] [185:214] [215:244]};

rd_diff_mean = NaN(6,N_LAT,N_LON);
for c=1:length(YRS_IND)
    rd_diff_mean(c,:,:) = nanmean(rd_diff(YRS_IND{c},:,:),1);
end

% Linear trend for plots; boxplots for periods.
rd_r45_trend = NaN(N_LAT,N_LON);
rd_r45_sig = NaN(N_LAT,N_LAT);
rd_r85_trend = NaN(N_LAT,N_LON);
rd_r85_sig = NaN(N_LAT,N_LON);
for i=1:N_LON
    for j=1:N_LAT
        [rd_r45_trend(j,i),rd_r45_sig(j,i)] = trenddata(rd(61:150,j,i));
        [rd_r85_trend(j,i),rd_r85_sig(j,i)] = trenddata(rd(155:244,j,i));
    end
end

data = rd_hst_mean;
map_title = 'Relative Difference Between LSF and GSI 1950-2005';
cb_type = 'div';
cb_color = 'RdBu';
cb_units = 'Difference in Days';
cb_flip = 'No Flip';
plotfscomponents(data,-30,30,6,map_title,cb_type,cb_color,cb_units,cb_flip)


data = rd_r85_trend;
map_title = 'Relative Difference Trend (2010-2099; RCP8.5)';
cb_type = 'div';
cb_color = 'RdBu';
cb_units = 'Trend over 2010-2099 in days';
cb_flip = 'Flip';
plotfscomponents(data*90,-40,40,10,map_title,cb_type,cb_color,cb_units,cb_flip)

% Calculate FSEI for 56 years 2044-2099.
YRS_IND = {[1:56] [95:150] [189:244]};

% Preallocate, iterate over periods.
fs_sum_v2 = NaN(3,N_LAT,N_LON);
fsei_v2 = NaN(3,N_LAT,N_LON);
for c=1:length(YRS_IND)
    fsei_v2(c,:,:) = (fs_sum(c,:,:) ./ sum(~isnan(lsf(YRS_IND{c},:,:)),1)) * 100;
end

% Plot future FSEI.
data = squeeze(fsei_v2(3,:,:));
map_title = 'False Spring Exposure Index (2044-2099; RCP8.5)';
cb_type = 'div';
cb_color = 'RdBu';
cb_units = 'Percent of Years';
cb_flip = 'Flip';
plotfscomponents(data,0,100,10,map_title,cb_type,cb_color,cb_units,cb_flip)

% Difference in trend, not trend in relative difference.
r45_trend_diff = lsf_trend - gsi_trend;
r85_trend_diff = lsf_r85_trend - gsi_r85_trend;

data = r45_trend_diff;
map_title = 'Difference in Trends (2010-2099; RCP4.5)';
cb_type = 'div';
cb_color = 'RdBu';
cb_units = 'Difference (days over 90 years)';
cb_flip = 'No Flip';
plotfscomponents(data*90,-60,60,15,map_title,cb_type,cb_color,cb_units,cb_flip)

% Make barplot figures showing percent of pixels on x-axis and FSEI on
% y-axis. Should be 2 rows, four columns corresponding to 2 RCPs and four
% time periods.
fsei_hst = squeeze(fsei(1,:,:));

test = numel(find(fsei>=50-1 & fsei<=50+1))  % Find all points between 49 and 50
grid_sum = numel(find(~isnan(fsei)))         % Sum all non-NaN grid points.
test_ans = (test/grid_sum) * 100             % Find percent.

for c=1:7
    t_fsei = squeeze(round(fsei(c,:,:)));
    t_sum = numel(find(~isnan(t_fsei)));
    for i=1:10
        n_fsei(c,i) = numel(find(t_fsei >= i-1 & t_fsei <= i+9));
        n_pxls(c,i) = (n_fsei(i) / t_sum) * 100;
    end
end

for c=1:7
    t_fsei = squeeze(round(fsei(c,:,:)));
    t_sum = numel(find(~isnan(t_fsei)));
    n_fsei(c,1) = numel(find(t_fsei>=0 & t_fsei<=9));
    n_fsei(c,2) = numel(find(t_fsei>=10 & t_fsei<=19));
    n_fsie(c,3) = numel(find(t_fsei>=20 & t_fsei<=29));
    n_fsei(c,4) = numel(find(t_fsei>=30 & t_fsei<=39));
    n_fsie(c,5) = numel(find(t_fsei>=40 & t_fsei<=49));
    n_fsei(c,6) = numel(find(t_fsei>=50 & t_fsei<=59));
    n_fsei(c,7) = numel(find(t_fsei>=60 & t_fsei<=69));
    n_fsei(c,8) = numel(find(t_fsei>=70 & t_fsei<=79));
    n_fsei(c,9) = numel(find(t_fsei>=80 & t_fsei<=89));
    n_fsei(c,10) = numel(find(t_fsei>=90 & t_fsei<=100));
    n_pxls(c,:) = (n_fsei(c,:) / t_sum) * 100;
    clear t_fsei t_sum
end


h = figure('Position',[100 100 1000 618]);
subaxis(2,3,1,'SpacingHoriz',0,'SpacingVert',0);
bar(n_pxls(2,:)','FaceColor',rgb('SteelBlue'),'LineStyle','None','BarWidth',1)
axis([0 10 1 50])
set(gca,'XTick',[1:1:9],...
        'XTickLabel','',...
        'YTick',[10:10:40],...
        'YTickLabel',(10:10:40))
ylabel({'Areal Percent'; 'RCP4.5'},'FontSize',14)
title({'2010-2039'},'FontSize',14)

subaxis(2,3,2,'SpacingHoriz',0,'SpacingVert',0);
bar(n_pxls(3,:)','FaceColor',rgb('SteelBlue'),'LineStyle','None','BarWidth',1)
axis([0 10 1 50])
set(gca,'XTick',[1:1:9],...
        'XTickLabel','',...
        'YTick',[10:10:40],...
        'YTickLabel','')
title({'Change in FSEI'; '2040-2069'},'FontSize',14)

subaxis(2,3,3,'SpacingHoriz',0,'SpacingVert',0);
bar(n_pxls(4,:)','FaceColor',rgb('SteelBlue'),'LineStyle','None','BarWidth',1)
axis([0 10 1 50])
set(gca,'XTick',[1:1:9],...
        'XTickLabel','',...
        'YTick',[10:10:40],...
        'YTickLabel','')
title({'2070-2099'},'FontSize',14)

subaxis(2,3,4,'SpacingHoriz',0,'SpacingVert',0);
bar(n_pxls(5,:)','FaceColor',rgb('IndianRed'),'LineStyle','None','BarWidth',1)
axis([0 10 1 50])
set(gca,'XTick',[1:1:9],...
        'XTickLabel',(10:10:90),...
        'YTick',[10:10:40],...
        'YTickLabel',(10:10:40))
ylabel({'Areal Percent'; 'RCP8.5'},'FontSize',14)

subaxis(2,3,5,'SpacingHoriz',0,'SpacingVert',0);
bar(n_pxls(6,:)','FaceColor',rgb('IndianRed'),'LineStyle','None','BarWidth',1)
axis([0 10 1 50])
set(gca,'XTick',[1:1:9],...
        'XTickLabel',(10:10:90),...
        'YTick',[10:10:40],...
        'YTickLabel','')
xlabel({'FSEI'},'FontSize',14)

subaxis(2,3,6,'SpacingHoriz',0,'SpacingVert',0);
bar(n_pxls(7,:)','FaceColor',rgb('IndianRed'),'LineStyle','None','BarWidth',1)
axis([0 10 1 50])
set(gca,'XTick',[1:1:9],...
        'XTickLabel',(10:10:90),...
        'YTick',[10:10:40],...
        'YTickLabel','')


h = figure('Position',[100 100 1000 618]);
subaxis(2,3,1,'SpacingHoriz',0,'SpacingVert',0);
bar(n_pxls(1,:)','FaceColor',rgb('DarkGray'),'LineStyle','None','BarWidth',1)
axis([0 10 1 50])
set(gca,'XTick',[1:1:9],...
        'XTickLabel',(10:10:90),...
        'YTick',[10:10:40],...
        'YTickLabel',(10:10:40))
ylabel({'Areal Percent'},'FontSize',14)
xlabel({'FSEI'},'FontSize',14)
title({'1950-2005'},'FontSize',14)



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
station_tmin = tmindata(:,103:158,:);
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





% Moving average test.
mdl_gsi_rshp = reshape(mdl_gsi,[30 585*1386]);

tic
for i=size(mdl_gsi_rshp,2)
    gsi_avg(:,i) = movingaverage(mdl_gsi(:,i),10);
end
toc

tic
for i=1:size(mdl_gsi_rshp,2)
    gsi_avg_v2(:,i) = movingaveragev2(mdl_gsi(:,i),10);
end
toc


output = tsmovavg(tsobj,'s',lag,dim)

tic
for i=1:size(mdl_gsi_rshp,2)
    gsi_avg_v2(:,i) = tsmovavg(mdl_gsi,'s',21,1);
end






















