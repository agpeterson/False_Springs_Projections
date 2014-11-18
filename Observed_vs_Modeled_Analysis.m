%%=============================================================================
% NAME:   Observed_vs_Modeled_Analysis.m
% AUTHOR: Alexander Peterson
% DATE:   21 Oct. 2014
% DESCR:  This script compares the observed GRIDMET data to the MACAv2-METDATA
%		  dataset.
% IN:     gridmet_19792009.mat; gsi.mat; lsf.mat
% OUT:    N/A
% CALLS:  
%==============================================================================

% Load datasets.
load gridmet_19792009
obs_gsi = double(gridmet_gsi);
obs_lsf = double(gridmet_lsf);
clear gridmet_gsi gridmet_lsf

% Load GSI and LSF using matfile pointers.
file = matfile('gsi.mat');
mdl_gsi = file.gsi_CONUS(30:59,:,:,6);
mdl_gsi = double(mdl_gsi);
clear file

file = matfile('lsf.mat');
mdl_lsf = file.lsf_CONUS(30:59,:,:,6);
mdl_lsf = double(mdl_lsf);
clear file

% Create constant for number of years, lat, lon.
N_LAT = 585;
N_LON = 1386;
N_YRS = 30;


%%=============================================================================
% Derive false springs for both datasets.
%==============================================================================
% Find years with false springs across all lat/lon points for MACA.
mdl_fs = NaN(N_YRS,N_LAT,N_LON);
for i=1:N_LON
    for j=1:N_LAT
        for k=1:N_YRS
            if isnan(mdl_lsf(k,j,i)) || isnan(mdl_gsi(k,j,i))
                mdl_fs(k,j,i) = NaN;
            elseif mdl_lsf(k,j,i) >= (7 + mdl_gsi(k,j,i))
                mdl_fs(k,j,i) = 1;
            else
                mdl_fs(k,j,i) = 0;
            end
        end
    end
end

% Find years with false springs across all lat/lon points for GRIDMET.
obs_fs = NaN(N_YRS,N_LAT,N_LON);
for i=1:N_LON
    for j=1:N_LAT
        for k=1:N_YRS
            if isnan(obs_lsf(k,j,i)) || isnan(obs_gsi(k,j,i))
                obs_fs(k,j,i) = NaN;
            elseif obs_lsf(k,j,i) >= (7 + obs_gsi(k,j,i))
                obs_fs(k,j,i) = 1;
            else
                obs_fs(k,j,i) = 0;
            end
        end
    end
end

% Derive FSEI for MACA.
mdl_fs_sum = squeeze(sum(mdl_fs(:,:,:),1));
mdl_lsf_sum = squeeze(sum(~isnan(mdl_lsf(:,:,:)),1));
mdl_fsei = (mdl_fs_sum(:,:) ./ 30) * 100;

obs_fs_sum = squeeze(sum(obs_fs(:,:,:),1));
obs_lsf_sum = squeeze(sum(~isnan(obs_lsf(:,:,:)),1));
obs_fsei = (obs_fs_sum(:,:) ./ 30) * 100;


%%=============================================================================
% Derive climatological normals for 1979-2009 and plot.
%==============================================================================
% Calculate mean.
obs_gsi_mean = squeeze(nanmean(obs_gsi,1));
obs_lsf_mean = squeeze(nanmean(obs_lsf,1));

mdl_gsi_mean = squeeze(nanmean(mdl_gsi,1));
mdl_lsf_mean = squeeze(nanmean(mdl_lsf,1));

% Flip GRIDMET columns to match MACA.
obs_gsi_mean = flipud(obs_gsi_mean);
obs_lsf_mean = flipud(obs_lsf_mean);
obs_fsei = flipud(obs_fsei);

% Reshape all to continuous series.
obs_gsi_mean_rshp = reshape(obs_gsi_mean,[1 N_LAT*N_LON]);
obs_lsf_mean_rshp = reshape(obs_lsf_mean,[1 N_LAT*N_LON]);
obs_fsei_rshp = reshape(obs_fsei,[1 N_LAT*N_LON]);

mdl_gsi_mean_rshp = reshape(mdl_gsi_mean,[1 N_LAT*N_LON]);
mdl_lsf_mean_rshp = reshape(mdl_lsf_mean,[1 N_LAT*N_LON]);
mdl_fsei_rshp = reshape(mdl_fsei,[1 N_LAT*N_LON]);


%%=============================================================================
% Calculate bivariate statistics.
%==============================================================================
% GSI
x1 = obs_gsi_mean_rshp;
y1 = mdl_gsi_mean_rshp;
[gsi_r gsi_p] = corrcoef(x1,y1,'rows','pairwise')
p1 = polyfit(x1(~isnan(x1)),y1(~isnan(y1)),1)
yfit1 = polyval(p1,x1);
tbl = table(obs_gsi_mean_rshp',mdl_gsi_mean_rshp',...
	'VariableNames',{'Obs','Mdl'});
mdl = fitlm(tbl)

% LSF
x2 = obs_lsf_mean_rshp;
y2 = mdl_lsf_mean_rshp;
[lsf_r lsf_p] = corrcoef(x2,y2,'rows','pairwise')
p2 = polyfit(x2(~isnan(y2)),y2(~isnan(y2)),1)
yfit2 = polyval(p2,x2);
tbl = table(obs_lsf_mean_rshp',mdl_lsf_mean_rshp',...
	'VariableNames',{'Obs','Mdl'});
mdl = fitlm(tbl)

% FSEI
x3 = obs_fsei_rshp;
y3 = mdl_fsei_rshp;
[fsei_r fsei_p] = corrcoef(x3,y3,'rows','pairwise')
p3 = polyfit(x3(~isnan(x3)),y3(~isnan(x3)),1)
yfit3 = polyval(p3,x3);
tbl = table(obs_fsei_rshp',mdl_fsei_rshp',...
	'VariableNames',{'Obs','Mdl'});
mdl = fitlm(tbl)


%%=============================================================================
% Plot bivariate analysis.
%==============================================================================
% GSI
subplot(3,1,1)
plot(x1,y1,'b+',x1,yfit1,'k','LineWidth',1.5)
title('Observed vs. Modeled GSI (1979-2009)')
xlabel('Observed DoY'); ylabel('Modeled DoY')
legend('Data','Regression','Location','Southeast')
grid on
axis([0 200 0 200])
set(gca,'YTick',[0:50:200])
text(5,150,{'R2: 0.997'; 'RMSE: 1.35'})

% LSF
subplot(3,1,2)
plot(x2,y2,'b+',x2,yfit2,'k','LineWidth',1.5)
title('Observed vs. Modeled LSF (1979-2009)')
xlabel('Observed DoY'); ylabel('Modeled DoY')
grid on
axis([0 200 0 200])
set(gca,'YTick',[0:50:200])
text(5,150,{'R2: 0.986'; 'RMSE: 4.1'})

% FSEI
subplot(3,1,3)
plot(x3,y3,'b+',x3,yfit3,'k','LineWidth',1.5)
title('Observed vs. Modeled FSEI (1979-2009)')
xlabel('Observed Percent'); ylabel('Modeled Percent')
grid on
text(3,75,{'R2: 0.719'; 'RMSE: 11.1'})


%%=============================================================================
% Derive differences and map.
%==============================================================================
% Take differences.
gsi_diff = mdl_gsi_mean - obs_gsi_mean;
lsf_diff = mdl_lsf_mean - obs_lsf_mean;
fsei_diff = mdl_fsei - obs_fsei;

% Load lat/lon.
load conus_grid
lon = lon - 360;

% GSI
data = gsi_diff;
map_title = '1979-2009 GSI Difference (Mdl - Obs)';
cb_type = 'div';
cb_color = 'RdBu';
cb_units = 'Number of Days';
cb_flip = 'No Flip';
mapgriddata(data,-7,7,2,lat,lon,map_title,...
    cb_type,cb_color,cb_units,cb_flip)

% LSF
data = lsf_diff;
map_title = '1979-2009 LSF Difference (Mdl - Obs)';
cb_type = 'div';
cb_color = 'RdBu';
cb_units = 'Number of Days';
cb_flip = 'No Flip';
mapgriddata(data,-14,14,4,lat,lon,map_title,...
    cb_type,cb_color,cb_units,cb_flip)

% FSEI
data = fsei_diff;
map_title = '1979-2009 FSEI Difference (Mdl - Obs)';
cb_type = 'div';
cb_color = 'RdBu';
cb_units = 'Percent';
cb_flip = 'No Flip';
mapgriddata(data,-50,50,10,lat,lon,map_title,...
    cb_type,cb_color,cb_units,cb_flip)


%%=============================================================================
% Plot distributions.
%==============================================================================
% GSI.
subaxis(3,2,1,'SH',0,'SV',0)
histx(obs_gsi_mean_rshp)
axis([0 200 0 15000])
grid on
set(gca,'XTickLabel',[],'YTick',[5000 10000])
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b','EdgeColor','b')
title('Observed (1979-2009)')
ylabel('GSI')

subaxis(3,2,2,'SH',0,'SV',0)
histx(mdl_gsi_mean_rshp)
axis([0 200 0 15000])
grid on
set(gca,'XTickLabel',[],'YTickLabel',[],'YTick',[5000 10000])
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b','EdgeColor','b')
title('Modeled (1979-2009)')

% LSF.
subaxis(3,2,3,'SH',0,'SV',0)
histx(obs_lsf_mean_rshp(obs_lsf_mean_rshp < 181))
axis([0 200 0 15000])
grid on
set(gca,'YTick',[5000 10000],'XTick',[50:50:150])
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b','EdgeColor','b')
ylabel('LSF'); xlabel('Day of Year')

subaxis(3,2,4,'SH',0,'SV',0)
histx(mdl_lsf_mean_rshp)
axis([0 200 0 15000])
grid on
set(gca,'YTick',[5000 10000],'XTick',[50:50:150],'YTickLabel',[])
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b','EdgeColor','b')
xlabel('Day of Year')

% FSEI.
subaxis(3,2,5,'SH',0,'SV',0,'PT',0.06)
histx(obs_fsei_rshp)
axis([0 100 0 40000])
grid on
set(gca,'YTick',[10000:10000:30000],'YTickLabel',[10000:10000:30000],'XTick',[25:25:75])
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b','EdgeColor','b')
ylabel('FSEI'); xlabel('Percent of Years')

subaxis(3,2,6,'SH',0,'SV',0,'PT',0.06)
histx(mdl_fsei_rshp)
axis([0 100 0 40000])
grid on
set(gca,'YTick',[10000:10000:30000],'YTickLabel',[],'XTick',[25:25:75])
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b','EdgeColor','b')
xlabel('Percent of Years')


%% maybe a ttest2?


for i=1:N_LON
    for j=1:N_LAT
        [H(j,i),P(j,i)] = ttest2(obs_gsi(:,j,i),mdl_gsi(:,j,i));
    end
end

load Data/conus_grid
lon = lon - 360;

% GSI
data = P;
map_title = 'TTest2 GSI';
cb_type = 'seq';
cb_color = 'Blues';
cb_units = 'Significance';
cb_flip = 'No Flip';
mapgriddata(data,0,1,.1,lat,lon,map_title,...
    cb_type,cb_color,cb_units,cb_flip)







%% New 
file = matfile('gridMET_gsi.mat');
gu = double(file.gsi_CONUS);

file = matfile('gridMET_vpd_gu.mat');
gu_vpd = double(file.gu_CONUS);

load('~/Dropbox/Workspace/ACSL/gridMET_Coords.mat')

% Climo.
gu_mean = squeeze(nanmean(gu,1));
gu_vpd_mean = squeeze(nanmean(gu_vpd,1));

% Map.
prj = 'Albers Equal-Area Conic';
map_title = 'Non-VPD Mean Green-up (1979-2012)'
cb_type = 'seq';
cb_color = 'Reds';
cb_units = 'Day of Year';
cb_flip = 'No Flip';
figure();
mapGriddedData(gu_mean,prj,1,180,20,lat,lon,2,2,map_title,cb_type,cb_color,...
    cb_units,cb_flip)


map_title = 'VPD Mean Green-up (1979-2012)';
mapGriddedData(gu_vpd_mean,prj,1,180,20,lat,lon,2,2,map_title,cb_type,cb_color,...
    cb_units,cb_flip)


% ttest
[h,p] = ttest2(gu,gu_vpd);
h = squeeze(h);
p = squeeze(p);

% Look at LSF.
file = matfile('gridMET_lsf.mat')
lsf = double(file.lsf_CONUS);

file = matfile('gridMET_vpd_lsf.mat')
lsf_vpd = double(file.lsf_CONUS);

% climo.
lsf_mean = squeeze(nanmean(lsf,1));
lsf_vpd_mean = squeeze(nanmean(lsf_vpd,1));


% Look at models - check for missing region.
mapGriddedData(lsf_mean,prj,1,180,20,lat,lon,2,2,map_title,cb_type,cb_color,...
    cb_units,cb_flip)


% models
PATH = '/media/alexander/Vault/gsi_'
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


m = matfile([PATH MDL_NAME{1}])
mdl1_gsi = m.gsi_CONUS;

lat2 = m.lat;
lon2 = m.lon;
mapGriddedData(squeeze(mdl1_gsi(50,:,:)),prj,1,180,20,lat2,lon2,2,2,map_title,cb_type,cb_color,...
    cb_units,cb_flip)





