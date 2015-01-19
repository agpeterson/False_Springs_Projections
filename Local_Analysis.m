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
VAR_NAME = {'gu_';...
			'lsf_'};
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
climo_ind = {hst_ind r45_fut1 r45_fut2 r85_fut1 r85_fut2};

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
addpath(genpath('../Function_Library'))

% Load MACA coords.
load MACA_Coords

% Load ecoregion masks.
load Ecoregion_Masks


%%=============================================================================
% Map individual model means and differences to check for errors.
%==============================================================================

% Iterate over number of models. Change _file variable.
for i=1:N_MDL
	data = squeeze(fsei_file.clm_mean(1,:,:,i));

	% Initialize variables.
	prj = 'Albers Equal-Area Conic';
	min_val = 0;
	max_val = 100;
	val_step = 10;
	lat_buffer = 2;
	lon_buffer = 2;
	cb_units = 'Percent of Years';
	cb_type = 'div';
	cb_color = 'RdBu';
	cb_flip = 'No Flip';

	% Change map_title and data variables then run.
	map_title = [MDL_NAME{i} ' FSEI 1950-2005'];
	figure('Position',[100 100 1000 618]);
	mapGriddedData(data,prj,min_val,max_val,val_step,...
	               lat,lon,lat_buffer,lon_buffer,...
	               map_title,cb_type,cb_color,cb_units,cb_flip)
end


%%=============================================================================
% Map multi-model mean climatologies.
%==============================================================================

% Iterate over number of models. Change _file variable and take mean over
% models.
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
	cb_type = 'div';
	cb_color = 'RdBu';
	cb_flip = 'No Flip';

	% Change map_title and data variables then run.
	map_title = ['FSEI 1950-2005'];
	figure('Position',[100 100 1000 618]);
	mapGriddedData(data,prj,min_val,max_val,val_step,...
	               lat,lon,lat_buffer,lon_buffer,...
	               map_title,cb_type,cb_color,cb_units,cb_flip)
end


%%=============================================================================
% Map multi-model future differences.
%==============================================================================

% Map climatological means.
clm = {'1950-2005' '2040-2069 RCP4.5' '2070-2099 RCP4.5' ...
	   '2040-2069 RCP8.5' '2070-2099 RCP8.5'};

% Iterate over number of models. Change _file variable and take mean over
% models.
for i=1:N_CLM

	% Load data.
	data = squeeze(gu_file.clm_diff(i,:,:,:));
	data = nanmean(double(data),3);

	% Initialize variables.
	prj = 'Albers Equal-Area Conic';
	min_val = -45;
	max_val = 0;
	val_step = 5;
	lat_buffer = 2;
	lon_buffer = 2;
	cb_units = 'Difference (Days)';
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
% Map multi-model future differences with model robustness measure.
% NOTE: FSEI needs to be tested using bootstrap. To be done on Thunder.
%==============================================================================

% Create hst_mean and fut_diff variables with climatologies for all models.
hst_mean = NaN(585,1386,20,3,'single');
hst_mean(:,:,:,1) = squeeze(gu_file.clm_mean(1,:,:,:));
hst_mean(:,:,:,2) = squeeze(lsf_file.clm_mean(1,:,:,:));
% hst_mean(:,:,:,3) = squeeze(fsei_file.clm_mean(1,:,:,:));

fut_diff = NaN(585,1386,20,3,'single');
fut_diff(:,:,:,1) = squeeze(gu_file.clm_diff(4,:,:,:));
fut_diff(:,:,:,2) = squeeze(lsf_file.clm_diff(4,:,:,:));
% fut_diff(:,:,:,3) = squeeze(fsei_file.clm_diff(4,:,:,:));

% Iterate over variables and call calcModelRobustness() function.
for i=1:3
	mdl_robustness(:,:,i) = calcModelRobustness(hst_mean(:,:,:,i),...
		fut_diff(:,:,:,i));
end

% Map climatological means.
var_name = {'GU' 'LSF' 'FSEI'};
for i=1:2

	% Load data.
	data = squeeze(nanmean(fut_diff(:,:,:,i),3));
	data(mdl_robustness(:,:,i)==0) = NaN;

	% Mapping variables.
	prj = 'Albers Equal-Area Conic';
	lat_buffer = 2;
	lon_buffer = 2;

	% Split on i (variable) where 1 and 2 share variables and 3 differs.
	if i < 3
		min_val = -45;
		max_val = 0;
		val_step = 5;
		cb_units = 'Difference (Days)';
		cb_type = 'seq';	%div
		cb_color = 'Reds';	%RdBu
		cb_flip = 'Flip';	% No Flip for FSEI
	else
		min_val = -50;
		max_val = 50;
		val_step = 10;
		cb_units = 'Difference (%)';
		cb_type = 'div';	%div
		cb_color = 'RdBu';	%RdBu
		cb_flip = 'No Flip';	% No Flip for FSEI
	end

	% Change map_title and data variables then run.
	map_title = ['Delta ' var_name{i} ' RCP8.5 2040-2069'];
	figure('Position',[100 100 1000 618]);
	mapGriddedData(data,prj,min_val,max_val,val_step,...
	               lat,lon,lat_buffer,lon_buffer,...
	               map_title,cb_type,cb_color,cb_units,cb_flip)

	% Clear data.
	clear data

end


%%=============================================================================
% Create ecoregions.
%==============================================================================

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

% Map ecoregion changes.
gu = squeeze(gu_file.clm_diff(4,:,:,:));
gu_mean = nanmean(double(gu),3);
conus_grid = NaN(585,1386,'single');
for i=1:19
	ecorgn_gu = gu_mean;
	ecorgn_gu(ecorgn_masks(:,:,i)==0) = NaN;
	ecorgn_mean = nanmean(reshape(ecorgn_gu,[585*1386 1]),1);
	conus_grid(ecorgn_masks(:,:,i)==1) = ecorgn_mean;
	clear ecorgn_gu ecorgn_mean
end

% Plot.
prj = 'Albers Equal-Area Conic';
min_val = -21;
max_val = 0;
val_step = 3;
lat_buffer = 2;
lon_buffer = 2;
cb_units = 'Difference (Days)';
cb_type = 'seq';	%div
cb_color = 'Reds';	%RdBu
cb_flip = 'Flip';	% No Flip for FSEI
map_title = 'Ecoregion Delta GU RCP8.5 2040-2069'
figure('Position',[100 100 1000 618]);
mapGriddedData(conus_grid,prj,min_val,max_val,val_step,...
	            lat,lon,lat_buffer,lon_buffer,...
	            map_title,cb_type,cb_color,cb_units,cb_flip)


%%=============================================================================
% Model spread boxplots.
%==============================================================================

% Create ecoregion means for all models.
for i=1:3

	% Switch cases on i to extract lat, lon, and mdl for variables.
    switch i,
        case 1, data = double(squeeze(gu_file.clm_mean(1,:,:,:)));
        case 2, data = double(squeeze(lsf_file.clm_mean(1,:,:,:)));
        case 3, data = double(squeeze(fsei_file.clm_mean(1,:,:,:)));
    end

	% Iterate over models and ecoregions, creating a temporary copy of model
	% data to mask over. Average over mask and hold.
	for j=1:20
		for k=1:19
			tmp = squeeze(data(:,:,j));
			tmp(ecorgn_masks(:,:,k)==0) = NaN;
			ecorgn_mean(k,j,i) = nanmean(reshape(tmp,[585*1386 1]),1);
		end
	end

end

% Set plotting variables.
symbols = {'o','o','o','o','o',...
		   's','s','s','s','s',...
		   'd','d','d','d','d',...
		   'p','p','p','p','p'};
mdl_colors = {'Purple','Red','Yellow','Green','Blue',...
			  'Purple','Red','Yellow','Green','Blue',...
			  'Purple','Red','Yellow','Green','Blue',...
			  'Purple','Red','Yellow','Green','Blue'};
ecorgn_names = {'Hot Cont.'...
				'Hot Cont. Mts.'...
				'Marine'...
				'Marine Mts.'...
				'Mediter.'...
				'Mediter. Mts.'...
				'Prairie'...
				'Savanna'...
				'Subtrop.'...
				'Subtrop. Mts.'...
				'Temper. Desert'...
				'Temper. Desert Mts.'...
				'Temper. Steppe'...
				'Temper. Steppe Mts.'...
				'Tropic. Subtrop. Desert'...
				'Tropic. Subtrop. Mts'...
				'Tropic. Subtrop. Steppe'...
				'Warm Cont.'...
				'Warm Cont. Mts.'};

% Create figure by calling subplot() and notBoxPlot3().
figure();

% GU
subaxis(1,3,1,'SpacingHoriz',0);
h = notBoxPlot3(ecorgn_mean(:,:,1)',[],0.3,'patch',symbols,mdl_colors);  % was 0.6
for k=1:19
    set(h(k).mu,'Color',rgb('Black'));
    set(h(k).semPtch,'FaceColor',rgb('DarkGray'));
    set(h(k).sdPtch,'FaceColor',rgb('DimGray'));
end
view(90,-90)	% Flip axes such that y is horizontal and x is vertical.
set(gca,'XDir','Reverse')	% Reverse the x-axis order.
set(gca,'YLim',[-50 50],...
		'YTick',(-40:20:40),...
		'YTickLabel',(-40:20:40),...
		'FontSize',10)
set(gca,'XTick',(1:19),'XTickLabel',ecorgn_names,'FontSize',10)
ylabel('Change (Days)','FontSize',12)
title('Delta GU 2040-2069')
grid off;
box on;
hold on;

% Model names.
for i=1:20
    plot(9.25+(i*.5),15,symbols{i},...
    	'MarkerFaceColor',rgb(mdl_colors{i}),...
        'MarkerEdgeColor',rgb('Black'));
    hold on
    text(9.25+(i*.5),17,MODEL{i},'fontsize',10);
end

% LSF
subaxis(1,3,2,'SpacingHoriz',0);
h = notBoxPlot3(ecorgn_mean(:,:,2)',[],0.3,'patch',symbols,mdl_colors);  % was 0.6
for k=1:19
    set(h(k).mu,'Color',rgb('Black'));
    set(h(k).semPtch,'FaceColor',rgb('DarkGray'));
    set(h(k).sdPtch,'FaceColor',rgb('DimGray'));
end
view(90,-90)	% Flip axes such that y is horizontal and x is vertical.
set(gca,'XDir','Reverse')	% Reverse the x-axis order.
set(gca,'YLim',[-50 50],...
		'YTick',(-40:20:40),...
		'YTickLabel',(-40:20:40),...
		'FontSize',10)
set(gca,'XTick',[1:19],'XTickLabel',[])
ylabel('Change (Days)','FontSize',12)
title('Delta LSF 2040-2069')
grid off;
box on;

% FSEI
subaxis(1,3,3,'SpacingHoriz',0);
h = notBoxPlot3(ecorgn_mean(:,:,3)',[],0.3,'patch',symbols,mdl_colors);  % was 0.6
for k=1:19
    set(h(k).mu,'Color',rgb('Black'));
    set(h(k).semPtch,'FaceColor',rgb('DarkGray'));
    set(h(k).sdPtch,'FaceColor',rgb('DimGray'));
end
view(90,-90)	% Flip axes such that y is horizontal and x is vertical.
set(gca,'XDir','Reverse')	% Reverse the x-axis order.
set(gca,'YLim',[-50 50],...
		'YTick',(-40:20:40),...
		'YTickLabel',(-40:20:40),...
		'FontSize',10)
set(gca,'XTick',[1:19],'XTickLabel',[])
ylabel('Change (%)','FontSize',12)
title('Delta FSEI 2040-2069')
grid off;
box on;

% Print statistics for each ecoregion.
for i=1:19
    mdl_max(i,:) = squeeze(max(ecorgn_mean(i,:,:)));
    mdl_mean(i,:) = squeeze(mean(ecorgn_mean(i,:,:)));
    mdl_min(i,:) = squeeze(min(ecorgn_mean(i,:,:)));
    mdl_stnd(i,:) = squeeze(std(ecorgn_mean(i,:,:)));
end
stats = [mdl_max,mdl_mean,mdl_min,mdl_stnd];

% Set variable names.
var_names = {'GU_max'...
             'LSF_max'...
             'FSEI_max'...
             'GU_mean'...
             'LSF_mean'...
             'FSEI_mean'...
             'GU_min'...
             'LSF_min'...
             'FSEI_min'...
             'GU_stnd'...
             'LSF_stnd'...
             'FSEI_stnd'};

% Create table and print.
summary_table = array2table(stats,'VariableNames',var_names,...
    'RowNames',ecorgn_names)
writetable(summary_table,'Hst_Stats','Delimiter',',','WriteRowNames',true)


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

% Look at histograms.
for i=1:20
    figure();
    histx(reshape(tmin_delta(:,:,i),[585*1386 1]))
end

% Flip to match GridMET coords.
tmin_delta = flipud(tmin_delta);

% Save to file to be uploaded to Thunder.
save('Tmin_Delta.mat','tmin_delta');

% Plot METDATA and sensitivity experiments.
file = matfile('/media/alexander/Vault/Bioclimate/GU_Sensitivity.mat')

% Plot.
prj = 'Albers Equal-Area Conic';
min_val = 1;
max_val = 181;
val_step = 20;
lat_buffer = 2;
lon_buffer = 2;
cb_units = 'Day of Year';
cb_type = 'div';	%div
cb_color = 'RdBu';	%RdBu
cb_flip = 'No Flip';	% No Flip for FSEI
map_title = 'GU Sensitivity';

for i=1:20
	disp(i)
	tmp = squeeze(nanmean(file.data(:,:,:,i),1));
	
	figure('Position',[100 100 1000 618]);
	mapGriddedData(tmp,prj,min_val,max_val,val_step,...
		            lat,lon,lat_buffer,lon_buffer,...
		            map_title,cb_type,cb_color,cb_units,cb_flip)
end





