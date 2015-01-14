%%=============================================================================
% NAME:   Model_Comparisons.m
% AUTHOR: Alexander Peterson
% DATE:   1 Dec. 2014
%
% DESC:   
% REF:	  None.
% NOTE:	  
%
% IN:     
% OUT:    
% CALL:   
%==============================================================================


%%=============================================================================
% Access LSF and GSI model output.
%==============================================================================

% Path suffix and prefix strings to be concatenated for model access.
PATH_PREFIX = '/media/alexander/Vault/Bioclimate/';	% Local
VARIABLE = {'gu_';...
			'lsf_'};
MODEL = {'bcc-csm1-1';...
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
FILE_EXT = '.mat';

% Create years variables.
YEARS = [1950:2005 2006:2099 2006:2099];

% Years to average over.
hst_start = 1950;
hst_end = 2005;

fut1_start = 2020;
fut1_end = 2059;

fut2_start = 2060;
fut2_end = 2099;

% Find start and end indices.
hst_start_ind = find(YEARS == hst_start);
hst_end_ind = find(YEARS == hst_end);

fut1_start_ind = find(YEARS == fut1_start);
fut1_end_ind = find(YEARS == fut1_end);

fut2_start_ind = find(YEARS == fut2_start);
fut2_end_ind = find(YEARS == fut2_end);

% Create RCP/year indices.
hst_ind = hst_start_ind:hst_end_ind;

r45_fut1 = fut1_start_ind(1):fut1_end_ind(1);
r45_fut2 = fut2_start_ind(1):fut2_end_ind(1);

r85_fut1 = fut1_start_ind(2):fut1_end_ind(2);
r85_fut2 = fut2_start_ind(2):fut2_end_ind(2);

% Create structure for climatology indices.
clm_ind = {hst_ind r45_fut1 r45_fut2 r85_fut1 r85_fut2};

% Create constants for number of years, lat, lon, variables, and models.
N_LAT = 585;
N_LON = 1386;
N_MDL = 20;
N_VAR = 2;
N_YRS = 244;
N_CLM = 5;		% Hst + R45 fut1 + R45 fut2 + R85 fut1 + R85 fut2

% Load data.
gu_file = matfile('/media/alexander/Vault/Bioclimate/GU.mat')
lsf_file = matfile('/media/alexander/Vault/Bioclimate/LSF.mat')
fsei_file = matfile('/media/alexander/Vault/Bioclimate/FSEI.mat')

% Add local paths.
addpath(genpath('../Data'))
addpath(genpath('../Plotting'))

% Load MACA coords.
load MACAv2METDATA_Coordinates
lon = lon - 360;


%%=============================================================================
% Map all models.
%==============================================================================

% Map historical means.
for i=1:N_MDL
	data = squeeze(fsei_file.clm_mean(5,:,:,i));

	% Initialize variables.
	prj = 'Albers Equal-Area Conic';
	min_val = 0;
	max_val = 100;
	val_step = 10;
	lat_buffer = 2;
	lon_buffer = 2;
	cb_units = 'Day of Year';
	cb_type = 'seq';
	cb_color = 'Blues';
	cb_flip = 'No Flip';

	% Change map_title and data variables then run.
	map_title = [MODEL{i} ' FSEI 1950-2005'];
	figure('Position',[100 100 1000 618]);
	mapGriddedData(data,prj,min_val,max_val,val_step,...
	               lat,lon,lat_buffer,lon_buffer,...
	               map_title,cb_type,cb_color,cb_units,cb_flip)
end


%%=============================================================================
% Map all models.
%==============================================================================

% Map historical means.
for i=1:N_CLM
	data = squeeze(fsei_file.clm_mean(i,:,:,:));
	data = nanmean(double(data),3);

	% Initialize variables.
	prj = 'Albers Equal-Area Conic';
	min_val = 0;
	max_val = 100;
	val_step = 10;
	lat_buffer = 2;
	lon_buffer = 2;
	cb_units = 'Percent of Years';
	cb_type = 'seq';
	cb_color = 'Blues';
	cb_flip = 'No Flip';

	% Change map_title and data variables then run.
	map_title = ['FSEI 1950-2005'];
	figure('Position',[100 100 1000 618]);
	mapGriddedData(data,prj,min_val,max_val,val_step,...
	               lat,lon,lat_buffer,lon_buffer,...
	               map_title,cb_type,cb_color,cb_units,cb_flip)
end


%%=============================================================================
% Create multi-model climatology maps.
%==============================================================================

% Map climatological means.
clm = {'1950-2005' '2040-2069 RCP4.5' '2070-2099 RCP4.5' ...
	   '2040-2069 RCP8.5' '2070-2099 RCP8.5'};
for i=1:N_CLM

	% Load data.
	data = squeeze(fsei_file.clm_diff(i,:,:,:));
	data = nanmean(double(data),3);

	% Initialize variables.
	prj = 'Albers Equal-Area Conic';
	min_val = -45;
	max_val = 0;
	val_step = 5;
	lat_buffer = 2;
	lon_buffer = 2;
	cb_units = 'Difference from Historic (Days)';
	cb_type = 'seq';	%div
	cb_color = 'Reds';	%RdBu
	cb_flip = 'Flip';	% No Flip for FSEI

	% Change map_title and data variables then run.
	map_title = ['Change in GU ' clm{i}];
	figure('Position',[100 100 1000 618]);
	mapGriddedData(data,prj,min_val,max_val,val_step,...
	               lat,lon,lat_buffer,lon_buffer,...
	               map_title,cb_type,cb_color,cb_units,cb_flip)
end


%%=============================================================================
% Create NCA regions and means.
%==============================================================================

% Create grid points.
[x,y] = meshgrid(lon,lat);
grid_points = [x(:) y(:)];

% Write to csv file.
dlmwrite('MACA_Coords.csv',grid_points,'precision','%.6f');

% Get all ecoregion text files in directory.
files = dir('../Data/Ecoregions_MACA/*.txt');
files = strvcat(files.name);

lon_grid = x;
lat_grid = y;

% Start local parallel pool.
gcp()

% Preallocate rgn_grid.
rgn_mask = NaN(585,1386,size(files,1));

% Call parallel function to strip files of whitespace, open file, and create
% index grid by calling findLatLonIndices() function (RENAME).
parfor i=1:size(files,1)
	file_name = strtrim(files(i,:));
	disp(file_name);
	rgn_coords = dlmread(file_name,',');
	rgn_mask(:,:,i) = createRegionMask(rgn_coords,lon_grid,lat_grid);
end

% test to see output.
load rgn_mask
tst = data;
tst(rgn_mask(:,:,1)==0) = NaN;

% Want to make blocks to be plotted, e.g., 1s for region 1, 2s for 2, etc...
eco_grid = NaN(585,1386,'single');
for i=1:19
	eco_grid(rgn_mask(:,:,i)==1) = i;
end
figure('Position',[100 100 1000 618]);
mapGriddedData(eco_grid,prj,1,19,1,lat,lon,lat_buffer,lon_buffer,...
	'CONUS Ecoregions','qual','Accent','Divison','No Flip')

% Print out one region historic GU.
data = squeeze(gu_file.clm_mean(1,:,:,:));
data = nanmean(double(data),3);
tst_grid = NaN(585,1386,'single');
for i=1:19
	tst = data;
	tst(rgn_mask(:,:,i)==0) = NaN;
	tst_mean = nanmean(reshape(tst,[585*1386 1]),1);
	tst_grid(rgn_mask(:,:,i)==1) = tst_mean;
	clear tst tst_mean
end
figure('Position',[100 100 1000 618]);
mapGriddedData(tst_grid,prj,1,181,20,lat,lon,lat_buffer,lon_buffer,...
	'Ecoregion Mean GU 1950-2005','seq','Blues','Day of Year','No Flip')


%%=============================================================================
% Model spread boxplots.
%==============================================================================

% Create ecoregion means for all models.
for i=1:3

	% Switch cases on i to extract lat, lon, and mdl for variables.
    switch i,
        case 1, data = double(squeeze(gu_file.clm_diff(5,:,:,:)));
        case 2, data = double(squeeze(lsf_file.clm_diff(5,:,:,:)));
        case 3, data = double(squeeze(fsei_file.clm_diff(5,:,:,:)));
    end

	% Iterate over models and ecoregions, creating a temporary copy of model
	% data to mask over. Mean over mask and store.
	for j=1:20

		for k=1:19
			tmp = squeeze(data(:,:,j));
			tmp(rgn_mask(:,:,k)==0) = NaN;
			rgn_mean(k,j,i) = nanmean(reshape(tmp,[585*1386 1]),1);
		end

	end

end

% Print statistics.
for i=1:19
	mdl_max(i,:) = squeeze(max(rgn_mean(i,:,:)));
	mdl_mean(i,:) = squeeze(mean(rgn_mean(i,:,:)));
	mdl_min(i,:) = squeeze(min(rgn_mean(i,:,:)));
end
stats = [mdl_max,mdl_mean,mdl_min];

% Write output to table.
row_names = {'HotContinentalDivision'...
			 'HotContinentalRegimeMountains'...
			 'MarineDivision'...
			 'MarineRegimeMountains'...
			 'MediterraneanDivision'...
			 'MediterraneanRegimeMountains'...
			 'PrairieDivision'...
			 'SavannaDivision'...
			 'SubtropicalDivision'...
			 'SubtropicalRegimeMountains'...
			 'TemperateDesertDivision'...
			 'TemperateDesertRegimeMountains'...
			 'TemperateSteppeDivision'...
			 'TemperateSteppeRegimeMountains'...
			 'TropicalSubtropicalDesertDivision'...
			 'TropicalSubtropicalRegimeMountains'...
			 'TropicalSubtropicalSteppeDivision'...
			 'WarmContinentalDivision'...
			 'WarmContinentalRegimeMountains'};
var_names = {'GU_max'...
			 'LSF_max'...
			 'FSEI_max'...
			 'GU_mean'...
			 'LSF_mean'...
			 'FSEI_mean'...
			 'GU_min'...
			 'LSF_min'...
			 'FSEI_min'};
summary_table = array2table(stats,'VariableNames',var_names,...
	'RowNames',row_names)
writetable(summary_table,'Fut2_R85_Stats','Delimiter',',','WriteRowNames',true)
	data = squeeze(fsei_file.clm_mean(5,:,:,i));


%%=============================================================================
% Sensitivity experiments.
%==============================================================================

% Look at delta tmin.
tmin_file = matfile('/media/alexander/Vault/Bioclimate/MACA_Tmin.mat')

% Calculate delta.
for i=1:20
	tmin_hst = tmin_file.data(:,:,1,i) - 273.15;
	tmin_fut = tmin_file.data(:,:,3,i) - 273.15;
	tmin_delta(:,:,i) = tmin_fut - tmin_hst;
end

for i=1:20
	figure();
	histx(reshape(tmin_delta(:,:,i),[585*1386 1]))
end

% Calculate multi-model mean tmin delta.
tmin_mdl_delta = squeeze(nanmean(tmin_delta,3));

% Load data.
gu_sen_file = matfile('/media/alexander/Vault/Bioclimate/GU_Sensitivity.mat')
gu_obs_file = matfile('/media/alexander/Vault/Bioclimate/GU_GridMET.mat')

lsf_sen_file = matfile('/media/alexander/Vault/Bioclimate/LSF_Sensitivity.mat')
lsf_obs_file = matfile('/media/alexander/Vault/Bioclimate/LSF_GridMET.mat')

% Take difference between climatological observed and sensitivity values.
gu_sen_mean = squeeze(nanmean(gu_sen_file.gu_CONUS,1));
gu_obs_mean = squeeze(nanmean(gu_obs_file.gu_CONUS,1));
gu_sen_delta = gu_sen_mean - gu_obs_mean;

lsf_sen_mean = squeeze(nanmean(lsf_sen_file.lsf_CONUS,1));
lsf_obs_mean = squeeze(nanmean(lsf_obs_file.lsf_CONUS,1));
lsf_sen_delta = lsf_sen_mean - lsf_obs_mean;

% Flip GridMET lats.
gu_sen_delta = flipud(gu_sen_delta);
lsf_sen_delta = flipud(lsf_sen_delta);

% Plot differences.
prj = 'Albers Equal-Area Conic';
min_val = -45;
max_val = 0;
val_step = 5;
lat_buffer = 2;
lon_buffer = 2;
cb_units = 'Difference (Days)';
cb_type = 'seq';
cb_color = 'Reds';
cb_flip = 'Flip';

map_title = 'Delta Mean LSF Date 1979-2012 (Sens.Exp.-Obs.)';
data = lsf_sen_delta;
figure('Position',[100 100 1000 618]);
mapGriddedData(data,prj,min_val,max_val,val_step,...
               lat,lon,lat_buffer,lon_buffer,...
               map_title,cb_type,cb_color,cb_units,cb_flip)

% Compare differences.
gu_fut_file = matfile('/media/alexander/Vault/Bioclimate/GU.mat')
lsf_fut_file = matfile('/media/alexander/Vault/Bioclimate/LSF.mat')

gu_clm_diff = squeeze(gu_fut_file.clm_diff(4,:,:,:));
gu_mdl_delta = squeeze(nanmean(gu_clm_diff,3));

lsf_clm_diff = squeeze(lsf_fut_file.clm_diff(4,:,:,:));
lsf_mdl_delta = squeeze(nanmean(lsf_clm_diff,3));

% Calculate correlations.
[r,p] = corrcoef(gu_sen_delta,gu_mdl_delta,'rows','pairwise')
[r,p] = corrcoef(lsf_sen_delta,lsf_mdl_delta,'rows','pairwise')

% Plot.
subplot(1,2,1)
scatter(reshape(gu_sen_delta,[585*1386 1]),reshape(gu_mdl_delta,[585*1386 1]),'+')
axis([-90 10 -90 10])
hold all
h = refline(1,0)
set(h,'Color','k','LineWidth',1.5)
set(gca,'XTick',[-80:20:0],...
	'XTickLabel',[-80:20:0],...
	'YTick',[-80:20:0],...
	'YTickLabel',[-80:20:0])
title('Delta GU RCP8.5 vs Sensitivity')
xlabel('Difference b/w Sens. & Obs.')
ylabel('Difference b/w Fut. & Hst.')
box on; grid on
text(-30,-70,{'Pearsons r: 0.9496'; 'p-value: 0'})
legend('Data','1:1 Reference Line','Location','Northwest')

subplot(1,2,2)
scatter(reshape(lsf_sen_delta,[585*1386 1]),reshape(lsf_mdl_delta,[585*1386 1]),'+')
axis([-90 10 -90 10])
hold all
h = refline(1,0)
set(h,'Color','k','LineWidth',1.5)
set(gca,'XTick',[-80:20:0],...
	'XTickLabel',[-80:20:0],...
	'YTick',[-80:20:0],...
	'YTickLabel',[-80:20:0])
title('Delta LSF RCP8.5 vs Sensitivity')
xlabel('Difference b/w Sens. & Obs.')
ylabel('Difference b/w Fut. & Hst.')
box on; grid on
text(-30,-70,{'Pearsons r: 0.8662'; 'p-value: 0'})

% Differences in the differences.
gu_deltas = gu_mdl_delta - gu_sen_delta;
lsf_deltas = lsf_mdl_delta - lsf_sen_delta;

% Plot differences.
prj = 'Albers Equal-Area Conic';
min_val = -15;
max_val = 15;
val_step = 3;
lat_buffer = 2;
lon_buffer = 2;
cb_units = 'Difference (Days)';
cb_type = 'div';
cb_color = 'RdBu';
cb_flip = 'No Flip';

map_title = 'Differences in GU Deltas (MACA - GridMET)';
data = gu_deltas;
figure('Position',[100 100 1000 618]);
mapGriddedData(data,prj,min_val,max_val,val_step,...
               lat,lon,lat_buffer,lon_buffer,...
               map_title,cb_type,cb_color,cb_units,cb_flip)




