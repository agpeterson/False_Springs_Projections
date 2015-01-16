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

% Load MACA grid.
load MACA_Coords

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

% Plot.
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

% Load ecoregion masks.
load Ecoregion_Masks

% Create ecoregion means for all models.
for i=1:3

	% Switch cases on i to extract lat, lon, and mdl for variables.
    switch i,
        case 1, data = double(squeeze(gu_file.clm_diff(4,:,:,:)));
        case 2, data = double(squeeze(lsf_file.clm_diff(4,:,:,:)));
        case 3, data = double(squeeze(fsei_file.clm_diff(4,:,:,:)));
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

