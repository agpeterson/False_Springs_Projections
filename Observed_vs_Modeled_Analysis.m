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
plot(x1,y1,'b+',x1,yfit1,'k','LineWidth',1.5)
title('Observed vs. Modeled GSI (1979-2009)')
xlabel('Observed DoY'); ylabel('Modeled DoY')
legend('Data','Regression','Location','Northwest')
grid on
axis([0 200 0 200])
set(gca,'YTick',[0:50:200])
text(5,150,{'R2: 0.997'; 'RMSE: 1.35'})

% LSF
plot(x2,y2,'b+',x2,yfit2,'k','LineWidth',1.5)
title('Observed vs. Modeled LSF (1979-2009)')
xlabel('Observed DoY'); ylabel('Modeled DoY')
legend('Data','Regression','Location','Northwest')
grid on
axis([0 200 0 200])
set(gca,'YTick',[0:50:200])
text(5,150,{'R2: 0.986'; 'RMSE: 4.1'})

% FSEI
plot(x3,y3,'b+',x3,yfit3,'k','LineWidth',1.5)
title('Observed vs. Modeled FSEI (1979-2009)')
xlabel('Observed Percent'); ylabel('Modeled Percent')
legend('Data','Regression','Location','Northwest')
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
mapgridded(data,-7,7,2,lat,lon,map_title,...
    cb_type,cb_color,cb_units,cb_flip)

% LSF
data = lsf_diff;
map_title = '1979-2009 LSF Difference (Mdl - Obs)';
cb_type = 'div';
cb_color = 'RdBu';
cb_units = 'Number of Days';
cb_flip = 'No Flip';
mapgridded(data,-14,14,4,lat,lon,map_title,...
    cb_type,cb_color,cb_units,cb_flip)

% FSEI
data = fsei_diff;
map_title = '1979-2009 FSEI Difference (Mdl - Obs)';
cb_type = 'div';
cb_color = 'RdBu';
cb_units = 'Percent';
cb_flip = 'No Flip';
mapgridded(data,-50,50,10,lat,lon,map_title,...
    cb_type,cb_color,cb_units,cb_flip)

