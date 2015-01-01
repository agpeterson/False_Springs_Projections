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
% Create historical climatology maps.
%==============================================================================

% Map historical means.
for i=1:N_MDL
	data = squeeze(gu_file.clm_mean(1,:,:,i));

	% Initialize variables.
	prj = 'Albers Equal-Area Conic';
	min_val = 1;
	max_val = 181;
	val_step = 20;
	lat_buffer = 2;
	lon_buffer = 2;
	cb_units = 'Day of Year';
	cb_type = 'seq';
	cb_color = 'Blues';
	cb_flip = 'No Flip';

	% Change map_title and data variables then run.
	map_title = [MODEL{i} ' Mean GU 1950-2005'];
	figure('Position',[100 100 1000 618]);
	mapGriddedData(data,prj,min_val,max_val,val_step,...
	               lat,lon,lat_buffer,lon_buffer,...
	               map_title,cb_type,cb_color,cb_units,cb_flip)
end


%%=============================================================================
% Create multi-model climatology maps.
%==============================================================================

% Map climatological means.
clm = {'1950-2005' '2020-2059 RCP4.5' '2060-2099 RCP4.5' ...
	   '2020-2059 RCP8.5' '2060-2099 RCP8.5'};
for i=1:N_CLM

	% Load data.
	data = squeeze(file.clm_diff(i,:,:,:));
	data = nanmean(double(data),3);

	% Initialize variables.
	prj = 'Albers Equal-Area Conic';
	min_val = -50;
	max_val = 50;
	val_step = 20;
	lat_buffer = 2;
	lon_buffer = 2;
	cb_units = 'Relative Difference (Percent)';
	cb_type = 'div';
	cb_color = 'RdBu';
	cb_flip = 'No Flip';

	% Change map_title and data variables then run.
	map_title = ['Difference in FSEI ' clm{i}];
	figure('Position',[100 100 1000 618]);
	mapGriddedData(data,prj,min_val,max_val,val_step,...
	               lat,lon,lat_buffer,lon_buffer,...
	               map_title,cb_type,cb_color,cb_units,cb_flip)
end


%%=============================================================================
% Create NCA regions and means.
%==============================================================================

% Write grid points to file.
[x,y] = meshgrid(lon,lat);
grid_points = [x(:) y(:)];
dlmwrite('MACA_Coords.csv',grid_points,'precision','%.6f');

% Load regional MACA coords.
nw_coords = dlmread('Northwest.txt',',', 2, 1);
sw_coords = dlmread('Southwest.txt',',', 2, 1);
gp_coords = dlmread('GreatPlains.txt',',', 2, 1);
mw_coords = dlmread('Midwest.txt',',', 2, 1);
ne_coords = dlmread('Northeast.txt',',', 2, 1);
se_coords = dlmread('Southeast.txt',',', 2, 1);

% Find regional coord indices.
[nw_lat,nw_lon] = findLatLonIndices(nw_coords,lat,lon);
[sw_lat,sw_lon] = findLatLonIndices(sw_coords,lat,lon);
[gp_lat,gp_lon] = findLatLonIndices(gp_coords,lat,lon);
[mw_lat,mw_lon] = findLatLonIndices(mw_coords,lat,lon);
[ne_lat,ne_lon] = findLatLonIndices(ne_coords,lat,lon);
[se_lat,se_lon] = findLatLonIndices(se_coords,lat,lon);

% Create cell array for regional averages.
region_lat = {nw_lat sw_lat gp_lat mw_lat ne_lat se_lat};
region_lon = {nw_lon sw_lon gp_lon mw_lon ne_lon se_lon};


%%=============================================================================
% Figure 1 boxplots - historical region means.
%==============================================================================

% Create regional and multi-model historical means.
for i=1:3

	% Switch cases on i to extract lat, lon, and mdl for variables.
	switch i,
		case 1, data = double(squeeze(gu_file.clm_mean(1,:,:,:)));
		case 2, data = double(squeeze(lsf_file.clm_mean(1,:,:,:)));
		case 3, data = double(squeeze(fsei_file.clm_mean(1,:,:,:)));
	end

	% Iterate over regions to derive regional means.
	for j=1:length(region_lat)

		% Index into data and extract region; reshape to array.
		data_rgn = data(region_lat{j},region_lon{j},:);

		% Set size variables.
		lat_size = size(data_rgn,1);
		lon_size = size(data_rgn,2);

		% Reshape vector to matrix.
		data_rshp = reshape(data_rgn,[lat_size*lon_size 20]);

		% Average over the region.
		rgn_mean(:,j,i) = nanmean(data_rshp,1)';
		
	end 	% k; 1:length(region_lat)

end 	% i; 1:3

% Multi-model mean.
rgn_mdl_mean = squeeze(nanmean(rgn_mean,1));

% Plot regional historic means for GU, LSF, and FSEI.
num_cats = 6;	% 6 regions.
num_subcats = 1;
x_labels = {'NW' 'SW' 'GP' 'MW' 'NE' 'SE'};
jitter = .5;
colors = {'DarkSlateGray' 'DarkGray' 'DarkSlateGray' 'DarkGray'};

for i=1:3
	
	% Create figures and squeeze data to be plotted.
	figure('Position',[100 100 1000 618])
	data = squeeze(rgn_mean(:,:,i));
	plotNotBoxPlot(data,num_cats,num_subcats,x_labels,jitter,colors)
	
	% Set y-axis limits.
	if i < 3
		set(gca,'YLim',[50 150],...
			'YTick',(70:20:130),...
			'YTickLabel',(70:20:130))
		legend('Mean','SEM','Model')
	else
		set(gca,'YLim',[0 100],...
			'YTick',(20:20:80),...
			'YTickLabel',(20:20:80))
		legend('Mean','SEM','Model')
	end

end

% Write output to table.
row_name = {'NW' 'SW' 'GP' 'MW' 'NE' 'SE'};
summary_table = array2table(rgn_mdl_mean,...
	'VariableNames',{'GU','LSF','FSEI'},'RowNames',row_name)
writetable(summary_table,'Hst_Means','Delimiter',',','WriteRowNames',true)


%%=============================================================================
% Figure 3 boxplots - regional future means.
%==============================================================================

% Preallocate variables.
rgn_mean = NaN(20,6,5,3);
rgn_mdl_mean = NaN(6,5,3);

% Create regional and multi-model means.
for i=1:3

	% Switch cases on i to extract lat, lon, and mdl for variables.
	switch i,
		case 1, data = double(squeeze(gu_file.clm_diff(:,:,:,:)));
		case 2, data = double(squeeze(lsf_file.clm_diff(:,:,:,:)));
		case 3, data = double(squeeze(fsei_file.clm_diff(:,:,:,:)));
	end

	% Iterate over climatologies.
	for j=1:N_CLM

		% Extract anomalies.
		data_clm = squeeze(data(j,:,:,:));

		% Iterate over regions to derive regional means.
		for k=1:length(region_lat)

			% Index into data and extract region; reshape to array.
			data_rgn = data_clm(region_lat{k},region_lon{k},:);

			% Set size variables.
			lat_size = size(data_rgn,1);
			lon_size = size(data_rgn,2);

			% Reshape vector to matrix.
			data_rshp = reshape(data_rgn,[lat_size*lon_size 20]);

			% Average over the region.
			rgn_mean(:,k,j,i) = nanmean(data_rshp,1)';
		
		end 	% k; 1:length(region_lat)

	end 	% j; 1:N_CLM
end 	% i; 1:3

% Multi-model mean.
rgn_mdl_mean = squeeze(nanmean(rgn_mean,1));

% Plot regional historic means for GU, LSF, and FSEI.
num_cats = 6;	% 6 regions.
num_subcats = 2;
x_labels = {'NW' 'SW' 'GP' 'MW' 'NE' 'SE'};
jitter = .6;
colors = {'DarkRed' 'IndianRed' 'DarkBlue' 'SteelBlue'};

% Create figures and squeeze data to be plotted.
figure('Position',[100 100 1000 618])
x = 1;
y = 2;
for i=1:3

	% Pull data; reshape to matrix, then rearrange rows.
	data = squeeze(rgn_mean(:,:,[2:5],i));

	data_1 = data(:,:,[1 3]);
	data_1 = reshape(data_1,[20 12]);
	data_1 = data_1(:,[1 7 2 8 3 9 4 10 5 11 6 12]);

	data_2 = data(:,:,[2 4]);
	data_2 = reshape(data_2,[20 12]);
	data_2 = data_2(:,[1 7 2 8 3 9 4 10 5 11 6 12]);

	% Plot
	subaxis(3,2,x,'SpacingHorizontal',0,'SpacingVertical',0);
	plotNotBoxPlot(data_1,num_cats,num_subcats,x_labels,jitter,colors)
	set(gca,'YLim',[-80 20],...
		'YTick',(-60:20:0),...
		'YTickLabel',(-60:20:0))

	if x == 1
		title('2020-2059')
		ylabel('Delta Days')
	elseif x == 3
		ylabel('Delta Days')
	elseif x == 5
		ylabel('Delta %')
	end		

	subaxis(3,2,y,'SpacingHorizontal',0,'SpacingVertical',0);
	plotNotBoxPlot(data_2,num_cats,num_subcats,x_labels,jitter,colors)
	set(gca,'YLim',[-80 20],...
		'YTick',(-60:20:0),...
		'YTickLabel',[])

	if y == 2
		title('2060-2099')
	end

	x = x + 2;
	y = y + 2;
	
end 	% i; 1:3

% Write output to table.
row_name = {'NW' 'SW' 'GP' 'MW' 'NE' 'SE'};
table_name = {'GU' 'LSF' 'FSEI'};
for i=1:3
	summary_table = array2table(squeeze(rgn_mdl_mean(:,2:5,i)),...
		'VariableNames',{'R45_Fut1','R45_Fut2','R85_Fut1','R85_Fut2'},...
		'RowNames',row_name)
	writetable(summary_table,['Fut_' table_name{i}],...
		'Delimiter',',','WriteRowNames',true)
end
