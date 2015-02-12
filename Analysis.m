%%============================================================================
%
% NAME:		Analysis.m
% AUTH:		Alexander Peterson
% DATE:		22 Jan. 2015
% DESC:		Script containing code analyzing GU, LSF, and FSEI.
%
%=============================================================================



%%============================================================================
% Initial variables and constants.
%-----------------------------------------------------------------------------

% Function and data paths.
% addpath(genpath('../Function_Library'))					% Local
% addpath(genpath('/media/alexander/Vault/Bioclimate/'))	% Local
addpath(genpath('False_Springs_Data/'))						% Thunder

% Variable and model names.
VAR_NAME = {'gu_'; 'lsf_'};
OBS_NAME = 'METDATA';
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


%-----------------------------------------------------------------------------
% Create constants for number of years, lat, lon, variables, and models.
%-----------------------------------------------------------------------------
N_LAT = 585;
N_LON = 1386;
N_MDL = 20;
N_VAR = 2;
N_YRS = 244;


%-----------------------------------------------------------------------------
% Matfile pointers and other files.
%-----------------------------------------------------------------------------

% Create pointer to model sensitivity matfiles.
sen_file = matfile('~/Model_Sensitivity.mat')
hst_file = matfile('~/Model_Historical.mat')
fut_file = matfile('~/Model_Future.mat')
obs_file = matfile('~/Observed.mat')

% Load extra files.
load('MACA_Coords.mat')
load('Ecoregion_Masks.mat')



%%============================================================================
% Create ecoregions.
%-----------------------------------------------------------------------------

% Create grid points.
[x,y] = meshgrid(lon,lat);
grid_points = [x(:) y(:)];

% Write to csv file.
dlmwrite('MACA_Coords.csv',grid_points,'precision','%.6f');

% Get all ecoregion text files in directory.
files = dir('../Data/Ecoregions_MACA/*.txt');
files = strvcat(files.name);

% Create lon/lat grids.
lon_grid = x;
lat_grid = y;

% Preallocate rgn_grid.
ecorgn_masks = NaN(585,1386,size(files,1));

% Strip files of whitespace, open file, and create index grid by calling 
% createRegionMask() function.
for i=1:size(files,1)
	file_name = strtrim(files(i,:));
	disp(file_name);
	ecorgn_coords = dlmread(file_name,',');
	ecorgn_masks(:,:,i) = createRegionMask(ecorgn_coords,lon_grid,lat_grid);
end
save('Ecoregion_Masks.mat','ecorgn_mask')

% Want to make blocks to be plotted, e.g., 1s for region 1, 2s for 2, etc...
ecorgn_grid = NaN(585,1386,'single');
for i=1:19
	ecorgn_grid(ecorgn_masks(:,:,i)==1) = i;
end

% Plot ecoregions.
prj = 'Albers Equal-Area Conic';
min_val = 1;
max_val = 20;
val_step = 1;
lat_buffer = 2;
lon_buffer = 2;
cb_units = 'Ecoregion Divison';
cb_type = 'qual';	%div
cb_color = 'Set1';	%RdBu
cb_flip = 'No Flip';	% No Flip for FSEI
map_title = 'CONUS Ecoregion Divisons'
figure('Position',[100 100 1000 618]);
mapGriddedData(ecorgn_grid,prj,min_val,max_val,val_step,...
	            lat,lon,lat_buffer,lon_buffer,...
	            map_title,cb_type,cb_color,cb_units,cb_flip)



%%============================================================================
% Create observed mean GU, LSF, and FSEI.
%-----------------------------------------------------------------------------

% Create matfile and preallocate GU, LSF, and FSEI values.
disp('Creating Observed_Historical.mat file...')
obs_file = matfile('Observed.mat','Writable',true);

% Concatenate to create filename for pointer.
disp('Accessing GU and LSF output from METDATA...')
gu_file = matfile('GU_METDATA.mat');
lsf_file = matfile('LSF_METDATA.mat');

% Load model data into memory.
disp('Loading data into memory...')
gu = double(gu_file.data);
lsf = double(lsf_file.data);

% Classify false springs and derive FSEI.
disp('Calling calcFSEI() function...')
[fs,fsei] = calcFSEI(gu,lsf,7);

% Calculate climatological mean GU and LSF.
disp('Taking climatological mean GU and LSF...')
gu_mean = squeeze(nanmean(gu,1));
lsf_mean = squeeze(nanmean(lsf,1));

% Write output to file.
disp('Writing output to file...')
obs_file.gu = single(gu);
obs_file.gu_mean = single(gu_mean);
obs_file.lsf = single(lsf);
obs_file.lsf_mean = single(lsf_mean);
obs_file.fs = single(fs);
obs_file.fsei = single(fsei);

% Clear memory.
clear gu_file lsf_file gu lsf fs fsei gu_mean lsf_mean



%%============================================================================
% Create modeled climatological mean GU, LSF, and FSEI.
%-----------------------------------------------------------------------------

% Create years constant, concatenating the historical 1950-2005 and two
% future 2006-2099 periods corresponding; the first future period corresponds
% to RCP4.5 and the second to RCP8.5.
YEARS = [1950:2005 2006:2099 2006:2099];

% Find year indices corresponding to climatology.
climatology = 2040:2069;
start_yrs_ind = find(YEARS == climatology(1));
end_yrs_ind = find(YEARS == climatology(end));
yrs_ind = [start_yrs_ind(2):end_yrs_ind(2)];

% Create matfile and preallocate GU, LSF, and FSEI values.
disp('Preallocating Model_Future.mat file...')
fut_file = matfile('Model_Future.mat','Writable',true);
fut_file.gu = NaN(30,N_LAT,N_LON,N_MDL,'single');
fut_file.gu_mean = NaN(N_LAT,N_LON,N_MDL,'single');
fut_file.lsf = NaN(30,N_LAT,N_LON,N_MDL,'single');
fut_file.lsf_mean = NaN(N_LAT,N_LON,N_MDL,'single');
fut_file.fs = NaN(30,N_LAT,N_LON,N_MDL,'single');
fut_file.fsei = NaN(N_LAT,N_LON,N_MDL,'single');

% Iterate over models for both variables, derive false springs, and take
% climatological means.
disp('Iterating over models to derive false springs and climatologies...')
for mdl=1:N_MDL

	tic

	% Concatenate to create filename for pointer.
	disp(['Accessing GU and LSF output from ' char(MDL_NAME{mdl}) '...'])
	gu_file = matfile([char(VAR_NAME{1}),char(MDL_NAME{mdl})]);
	lsf_file = matfile([char(VAR_NAME{2}),char(MDL_NAME{mdl})]);
	
	% Load model data into memory.
	disp('Loading data into memory...')
	gu = double(gu_file.gsi_CONUS(yrs_ind,:,:));
	lsf = double(lsf_file.lsf_CONUS(yrs_ind,:,:));

	% Classify false springs and derive FSEI.
	disp('Calling calcFSEI() function...')
	[fs,fsei] = calcFSEI(gu,lsf,7);

	% Calculate climatological mean GU and LSF.
	disp('Taking climatological mean GU and LSF...')
	gu_mean = squeeze(nanmean(gu,1));
	lsf_mean = squeeze(nanmean(lsf,1));

	% Write output to file.
	disp('Writing output to file...')
	fut_file.gu(:,:,:,mdl) = single(gu);
	fut_file.gu_mean(:,:,mdl) = single(gu_mean);
	fut_file.lsf(:,:,:,mdl) = single(lsf);
	fut_file.lsf_mean(:,:,mdl) = single(lsf_mean);
	fut_file.fs(:,:,:,mdl) = single(fs);
	fut_file.fsei(:,:,mdl) = single(fsei);

	% Clear memory.
	clear gu_file lsf_file gu lsf fs fsei gu_mean lsf_mean

	toc

end 	% mdl; 1:N_MDL



%%============================================================================
% Map historical and observed climatologies.
%-----------------------------------------------------------------------------

% Create multi-model mean historical values to be mapped and plotted.
gu_mdl = squeeze(nanmean(double(hst_file.gu_mean),3));
lsf_mdl = squeeze(nanmean(double(hst_file.lsf_mean),3));
fsei_mdl = squeeze(nanmean(double(hst_file.fsei),3));

% Load observed values into memory.
gu_obs = double(flipud(obs_file.gu_mean));
lsf_obs = double(flipud(obs_file.lsf_mean));
fsei_obs = double(flipud(obs_file.fsei));


%-----------------------------------------------------------------------------
% Map.
%-----------------------------------------------------------------------------

% Set mapGriddedData() function variables.
prj = 'Albers Equal-Area Conic';
min_val = 1;
lat_buffer = 2;
lon_buffer = 2;
cb_flip = 'No Flip';

% GU and LSF variables.
max_val = 181;
val_step = 20;
cb_units = 'Day of Year';
cb_type = 'div';
cb_color = 'RdBu';

% FSEI variables.
max_val = 100;
val_step = 10;
cb_units = 'Percent of Years';
cb_type = 'seq';
cb_color = 'Blues';

% GU.
map_title = ['Multi-Model Mean GU Date 1950-2005'];
figure('Position',[100 100 1000 618]);
mapGriddedData(gu_mdl,prj,min_val,max_val,val_step,...
               lat,lon,lat_buffer,lon_buffer,...
               map_title,cb_type,cb_color,cb_units,cb_flip)

% LSF.
map_title = ['Multi-Model Mean LSF Date 1950-2005'];
figure('Position',[100 100 1000 618]);
mapGriddedData(lsf_mdl,prj,min_val,max_val,val_step,...
               lat,lon,lat_buffer,lon_buffer,...
               map_title,cb_type,cb_color,cb_units,cb_flip)

% FSEI.
map_title = ['Multi-Model Mean FSEI 1950-2005'];
figure('Position',[100 100 1000 618]);
mapGriddedData(fsei_mdl,prj,min_val,max_val,val_step,...
               lat,lon,lat_buffer,lon_buffer,...
               map_title,cb_type,cb_color,cb_units,cb_flip)



%%============================================================================
% Compare modeled historical 1950-2005 to observed 1980-2009 values.
% Create linear models for r2 and RMSE to be plotted. First, reshape data
% into vector, then run correlations, fit, and evaluate. Finally, create table
% then linear model.
%-----------------------------------------------------------------------------

% GU
x = reshape(gu_mdl,[585*1386 1]);
y = reshape(gu_obs,[585*1386 1]);

% LSF
x = reshape(lsf_mdl,[585*1386 1]);
y = reshape(lsf_obs,[585*1386 1]);

% FSEI
x = reshape(fsei_mdl,[585*1386 1]);
y = reshape(fsei_obs,[585*1386 1]);

% Statistical body.
[r,p] = corrcoef(x,y,'rows','pairwise')
pfit = polyfit(x(~isnan(x)),y(~isnan(y)),1)
yfit = polyval(pfit,y);
tbl = table(x,y,'VariableNames',{'Mdl','Obs'});

% Create linear model.
gu_lm = fitlm(tbl)
lsf_lm = fitlm(tbl)
fsei_lm = fitlm(tbl)

% Clear variables.
clear r p pfit yfit tbl


%-----------------------------------------------------------------------------
% Plot linear models.
%-----------------------------------------------------------------------------

% GU
figure('Position',[100 100 618 618]);
plot(x,y,'+')
axis([0 200 0 200])
hold all
h = refline(1,0)
set(h,'Color','k','LineWidth',1.5)
set(gca,'XTick',[20:40:180],...
    'XTickLabel',[20:40:180],...
    'YTick',[20:40:180],...
    'YTickLabel',[20:40:180])
title({'Modeled 1950-2005 vs Observed 1980-2009'; 'Mean GU'},'FontSize',18)
xlabel('Modeled DoY','FontSize',16)
ylabel('Observed DoY','FontSize',16)
box on; grid on
text(135,40,{'Adjusted r^2: 0.997'; 'RMSE: 1.45'},'FontSize',14)
legend('Data','1:1 Reference Line','Location','Northwest')
legend('boxoff')

% LSF
figure('Position',[100 100 618 618]);
plot(x,y,'+')
axis([0 200 0 200])
hold all
h = refline(1,0)
set(h,'Color','k','LineWidth',1.5)
set(gca,'XTick',[20:40:180],...
    'XTickLabel',[20:40:180],...
    'YTick',[20:40:180],...
    'YTickLabel',[20:40:180])
title({''; 'Mean LSF'},'FontSize',18)
xlabel('Modeled DoY','FontSize',16)
ylabel('Observed DoY','FontSize',16)
box on; grid on
text(135,40,{'Adjusted r^2: 0.997'; 'RMSE: 2.13'},'FontSize',14)

% FSEI
figure('Position',[100 100 618 618]);
plot(x,y,'+')
axis([0 100 0 100])
hold all
h = refline(1,0)
set(h,'Color','k','LineWidth',1.5)
set(gca,'XTick',[10:20:90],...
    'XTickLabel',[10:20:90],...
    'YTick',[10:20:90],...
    'YTickLabel',[10:20:90])
title('FSEI','FontSize',18)
xlabel('Modeled %','FontSize',16)
ylabel('Observed %','FontSize',16)
box on; grid on
text(65,20,{'Adjusted r^2: 0.871'; 'RMSE: 7.37'},'FontSize',14)



%%============================================================================
% Model bias from observed.
%-----------------------------------------------------------------------------

% Calculate model bias by taking the difference between modeled and observed.
for mdl=1:N_MDL
	gu_hst_bias(:,:,mdl) = double(hst_file.gu_mean(:,:,mdl)) - gu_obs;
	lsf_hst_bias(:,:,mdl) = double(hst_file.lsf_mean(:,:,mdl)) - lsf_obs;
	fsei_hst_bias(:,:,mdl) = double(hst_file.fsei(:,:,mdl)) - fsei_obs;
end

% Calculate multi-model mean (mmm) bias.
gu_mmm_hst_bias = squeeze(nanmean(gu_hst_bias,3));
lsf_mmm_hst_bias = squeeze(nanmean(lsf_hst_bias,3));
fsei_mmm_hst_bias = squeeze(nanmean(fsei_hst_bias,3));

% Save biases to matfile.
gu_hst_bias = single(gu_hst_bias);
gu_mmm_hst_bias = single(gu_mmm_hst_bias);
lsf_hst_bias = single(lsf_hst_bias);
lsf_mmm_hst_bias = single(lsf_mmm_hst_bias);
fsei_hst_bias = single(fsei_hst_bias);
fsei_mmm_hst_bias = single(fsei_mmm_hst_bias);
save('Model_Historical_Biases.mat','gu_hst_bias','gu_mmm_hst_bias',...
	 'lsf_hst_bias','lsf_mmm_hst_bias','fsei_hst_bias','fsei_mmm_hst_bias')


%-----------------------------------------------------------------------------
% Run bootstrap method for difference/significance.
%-----------------------------------------------------------------------------

% Create matfile pointer to historical bias file.
hst_bias_file = matfile('Model_Historical_Biases.mat','Writable',true)
hst_bias_file.bootstrp = NaN(N_LAT,N_LON,N_MDL,'single');
hst_bias_file.signif = NaN(N_LAT,N_LON,N_MDL,'single');

% Load historical FS variable into memory.
fs_hst = double(hst_file.fs);
fs_obs = double(obs_file.fs);

% Permute and flip fs_obs to correct dimensions.
fs_obs = fs_obs(:,end:-1:1,:);

% Set yr_ind to random 30 years.
yr_ind = randi([1 56],1,30);

% Iterate over all models to call bootstrap function.
for mdl=1:N_MDL

	disp({'Model: ', mdl})

	% Pull false springs per model with random 30 years.
	disp('Accessing model false springs data...')
	fs_tmp = squeeze(fs_hst(yr_ind,:,:,mdl));

	% Call function.
	disp('Calling calcFSEIChangeSignif() function...')
	[fsei_bootstrp,fsei_signif] = calcFSEIChangeSignif(fs_tmp,fs_obs);

	% Write change and significance to file.
	disp('Writing output to file...')
	hst_bias_file.bootstrp(:,:,mdl) = single(fsei_bootstrp);
	hst_bias_file.signif(:,:,mdl) = single(fsei_signif);

	% Clear memory.
	clear fs_tmp fsei_bootstrp fsei_signif

end


%-----------------------------------------------------------------------------
% Map multi-model mean biases.
%-----------------------------------------------------------------------------

% Set mapGriddedData() function variables.
prj = 'Albers Equal-Area Conic';
lat_buffer = 2;
lon_buffer = 2;
cb_flip = 'No Flip';

% GU and LSF variables.
min_val = -10;
max_val = 10;
val_step = 4;
cb_units = 'Bias (Days)';
cb_type = 'div';
cb_color = 'RdBu';

% FSEI variables.
min_val = -28;
max_val = 28;
val_step = 8;
cb_units = 'Bias (Relative Percent)';
cb_type = 'div';
cb_color = 'RdBu';

% GU.
map_title = ['GU Bias (Modeled Historical - Observed)'];
figure('Position',[100 100 1000 618]);
mapGriddedData(double(gu_mmm_hst_bias),prj,min_val,max_val,val_step,...
               lat,lon,lat_buffer,lon_buffer,...
               map_title,cb_type,cb_color,cb_units,cb_flip)

% LSF.
map_title = ['LSF Bias (Modeled Historical - Observed)'];
figure('Position',[100 100 1000 618]);
mapGriddedData(double(lsf_mmm_hst_bias),prj,min_val,max_val,val_step,...
               lat,lon,lat_buffer,lon_buffer,...
               map_title,cb_type,cb_color,cb_units,cb_flip)

% FSEI.
map_title = ['FSEI Bias (Modeled Historical - Observed)'];
figure('Position',[100 100 1000 618]);
mapGriddedData(double(fsei_mmm_hst_bias),prj,min_val,max_val,val_step,...
               lat,lon,lat_buffer,lon_buffer,...
               map_title,cb_type,cb_color,cb_units,cb_flip)


%-----------------------------------------------------------------------------
% Explore FSEI bias using bootstrap sampling.
%-----------------------------------------------------------------------------

% Multi-model mean change.
mmm_btstrp = squeeze(nanmean(double(bootstrp),3)) * 100;

% Plot relative frequency.
figure('Position',[100 100 618 618]);
bins = -100:1:100;
[count,bin] = hist(reshape(mmm_btstrp,[N_LAT*N_LON 1]),bins); 
h = bar(bin,count/numel(mmm_btstrp),'k');
xlim([-100 100])
ylim([0 .2])
grid on
title('FSEI Bias (Bootstrap Resampling)')
xlabel('Bias (Relative Percent)'); ylabel('Relative Frequency')

% Set mapGriddedData() function variables.
prj = 'Albers Equal-Area Conic';
lat_buffer = 2;
lon_buffer = 2;
cb_flip = 'No Flip';
min_val = -28;
max_val = 28;
val_step = 8;
cb_units = 'Bias (Relative Percent)';
cb_type = 'div';
cb_color = 'RdBu';

% FSEI.
map_title = ['FSEI Bias (Bootstrap Resampling)'];
figure('Position',[100 100 1000 618]);
mapGriddedData(mmm_btstrp,prj,min_val,max_val,val_step,...
               lat,lon,lat_buffer,lon_buffer,...
               map_title,cb_type,cb_color,cb_units,cb_flip)


%-----------------------------------------------------------------------------
% Plot bias relative frequency distributions.
%-----------------------------------------------------------------------------

% GU.
figure('Position',[100 100 618 618]);
bins = -20:1:20;
[count,bin] = hist(reshape(gu_mmm_hst_bias,[N_LAT*N_LON 1]),bins);
h = bar(bin,count/numel(gu_mmm_hst_bias),'k');
xlim([-20 20])
ylim([0 .25])
grid on
title('GU Bias')
xlabel('Bias (Days)'); ylabel('Relative Frequency')

% LSF.
figure('Position',[100 100 618 618]);
bins = -20:1:20;
[count,bin] = hist(reshape(lsf_mmm_hst_bias,[N_LAT*N_LON 1]),bins); 
h = bar(bin,count/numel(lsf_mmm_hst_bias),'k');
xlim([-20 20])
ylim([0 .25])
grid on
title('LSF Bias')
xlabel('Bias (Days)'); ylabel('Relative Frequency')

% FSEI.
figure('Position',[100 100 618 618]);
bins = -40:1:40;
[count,bin] = hist(reshape(fsei_mmm_hst_bias,[N_LAT*N_LON 1]),bins); 
h = bar(bin,count/numel(fsei_mmm_hst_bias),'k');
xlim([-40 40])
ylim([0 .2])
grid on
title('FSEI Bias (Relative Difference)')
xlabel('Bias (Relative Percent)'); ylabel('Relative Frequency')



%%============================================================================
% Derive future change from 2040-2069 climatologies for GU and LSF.
%-----------------------------------------------------------------------------

% Load data into memory.
gu_hst = double(hst_file.gu_mean);
gu_fut = double(fut_file.gu_mean);

lsf_hst = double(hst_file.lsf_mean);
lsf_fut = double(fut_file.lsf_mean);

% Take difference for each model between future and historical values.
gu_mdl_delta = gu_fut - gu_hst;
lsf_mdl_delta = lsf_fut - lsf_hst;

% Create multi-model mean deltas.
gu_mmm_delta = nanmean(gu_mdl_delta,3);
lsf_mmm_delta = nanmean(lsf_mdl_delta,3);


%-----------------------------------------------------------------------------
% Map multi-model mean deltas.
%-----------------------------------------------------------------------------

% Set mapGriddedData() function variables.
prj = 'Albers Equal-Area Conic';
lat_buffer = 2;
lon_buffer = 2;

% GU and LSF variables.
min_val = 0;
max_val = 45;
val_step = 5;
cb_units = 'Advancement (Days)';
cb_type = 'seq';
cb_color = 'Reds';
cb_flip = 'No Flip';

% GU.
map_title = ['Change in GU Date RCP8.5 2040-2069'];
figure('Position',[100 100 1000 618]);
mapGriddedData(abs(gu_mmm_delta),prj,min_val,max_val,val_step,...
               lat,lon,lat_buffer,lon_buffer,...
               map_title,cb_type,cb_color,cb_units,cb_flip)

% LSF.
map_title = ['Change in LSF Date RCP8.5 2040-2069'];
figure('Position',[100 100 1000 618]);
mapGriddedData(abs(lsf_mmm_delta),prj,min_val,max_val,val_step,...
               lat,lon,lat_buffer,lon_buffer,...
               map_title,cb_type,cb_color,cb_units,cb_flip)



%%============================================================================
% Sensitivity analysis - compare model historical vs. experiments.
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% Create tmin variance deltas.
%-----------------------------------------------------------------------------

% Look at delta tmin.
file = matfile('/media/alexander/Vault/Bioclimate/MACA_Variance.mat')

% Calculate delta.
for i=1:N_MDL
    tmin_stnd_delta(:,:,i) = (fut_stnd(:,:,i) - hst_stnd(:,:,i)) ./ hst_stnd(:,:,i);
end

% Look at histograms.
for i=1:N_MDL
    figure();
    histx(reshape(tmin_stnd_delta(:,:,i),[585*1386 1]))
end

% Save to file to be uploaded to Thunder.
save('Tmin_Variance_Delta.mat','tmin_stnd_delta');
















% Load variables into memory.
gu_exp = double(exp_file.gu_mean);
lsf_exp = double(exp_file.lsf_mean);

% Find difference between modeled historical and sensitivty values.
for i=1:N_MDL
	gu_exp_delta(:,:,i) = gu_exp(:,:,i) - gu_hst(:,:,i);
	lsf_exp_delta(:,:,i) = lsf_exp(:,:,i) - lsf_hst(:,:,i);
end

% Take multi-model mean.
gu_mmm_exp_delta = squeeze(nanmean(gu_exp_delta,3));
lsf_mmm_exp_delta = squeeze(nanmean(lsf_exp_delta,3));


%-----------------------------------------------------------------------------
% Map multi-model mean deltas.
%-----------------------------------------------------------------------------

% Set mapGriddedData() function variables.
prj = 'Albers Equal-Area Conic';
lat_buffer = 2;
lon_buffer = 2;

% GU and LSF variables.
min_val = 0;
max_val = 45;
val_step = 5;
cb_units = 'Delta (Days)';
cb_type = 'seq';
cb_color = 'Reds';
cb_flip = 'No Flip';

% GU.
map_title = ['Change in GU Date Sensitivity 1950-2005'];
figure('Position',[100 100 1000 618]);
mapGriddedData(abs(gu_mmm_exp_delta),prj,min_val,max_val,val_step,...
               lat,lon,lat_buffer,lon_buffer,...
               map_title,cb_type,cb_color,cb_units,cb_flip)

% LSF.
map_title = ['Change in LSF Date Sensitivity 1950-2005'];
figure('Position',[100 100 1000 618]);
mapGriddedData(abs(lsf_mmm_exp_delta),prj,min_val,max_val,val_step,...
               lat,lon,lat_buffer,lon_buffer,...
               map_title,cb_type,cb_color,cb_units,cb_flip)


%-----------------------------------------------------------------------------
% Create linear models to compare deltas.
%-----------------------------------------------------------------------------

% GU deltas.
x1 = reshape(gu_mmm_delta,[585*1386 1]);
y1 = reshape(gu_mmm_exp_delta,[585*1386 1]);

% LSF
x2 = reshape(lsf_mmm_delta,[585*1386 1]);
y2 = reshape(lsf_mmm_exp_delta,[585*1386 1]);

% FSEI
x3 = reshape(fsei_mdl,[585*1386 1]);
y3 = reshape(fsei_obs,[585*1386 1]);

% Statistical body.
[r,p] = corrcoef(x2,y2,'rows','pairwise')
pfit = polyfit(x2(~isnan(x2)),y2(~isnan(y2)),1)
yfit = polyval(pfit,y2);
tbl = table(x2,y2,'VariableNames',{'Fut','Exp'});

% Create linear model.
gu_lm = fitlm(tbl)
lsf_lm = fitlm(tbl)
fsei_lm = fitlm(tbl)

% Clear variables.
clear r p pfit yfit tbl


%-----------------------------------------------------------------------------
% Plot linear models using color-coded scatterplot.
%-----------------------------------------------------------------------------

% Reshape ecoregions to continuous.
rgn = reshape(ecorgn_masks,[585*1386 19]);
rgn_colors = colormap(jet(19));

% GU deltas.
figure('Position',[100 100 618 618]);
plot(x1,y1,'+')
axis([-100 10 -100 10])
hold on
h = refline(1,0)
set(h,'Color','k','LineWidth',1.5)
set(gca,'XTick',[-80:20:0],...
    'XTickLabel',[-80:20:0],...
    'YTick',[-80:20:0],...
    'YTickLabel',[-80:20:0])
title({'Future vs Sensitivity Deltas'; 'GU'},'FontSize',16)
xlabel('Future Deltas','FontSize',14)
ylabel('Sensitivity Deltas','FontSize',14)
box on; grid on
text(-35,-70,{'Adjusted r2: 0.939'; 'RMSE: 1.66'},'FontSize',14)
legend('Data','1:1 Reference Line','Location','Northwest')
legend('boxoff')

% LSF
figure('Position',[100 100 618 618]);
plot(x2,y2,'+')
axis([-100 10 -100 10])
hold all
h = refline(1,0)
set(h,'Color','k','LineWidth',1.5)
set(gca,'XTick',[-80:20:0],...
    'XTickLabel',[-80:20:0],...
    'YTick',[-80:20:0],...
    'YTickLabel',[-80:20:0])
title({''; 'Mean LSF'},'FontSize',16)
xlabel('Future Deltas','FontSize',14)
ylabel('Sensitivity Deltas','FontSize',14)
box on; grid on
text(-35,-70,{'Adjusted r2: 0.88'; 'RMSE: 2.25'},'FontSize',14)

% FSEI
figure('Position',[100 100 618 618]);
plot(x3,y3,'+')
axis([0 100 0 100])
hold all
h = refline(1,0)
set(h,'Color','k','LineWidth',1.5)
set(gca,'XTick',[10:20:90],...
    'XTickLabel',[10:20:90],...
    'YTick',[10:20:90],...
    'YTickLabel',[10:20:90])
title('FSEI')
xlabel('Modeled %')
ylabel('Observed %')
box on; grid on
text(70,20,{'Adjusted r2: 0.871'; 'RMSE: 7.37'})


