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
% Create NCA regions.
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

% Plot test region (NW) using plotNotBoxPlot().
num_cats = 2;
num_subcats = 2;
x_labels = {'2020-2059' '2060-2099'};

data = squeeze(rgn_mean(:,1,2:end,1)); % All models, NW, futures, GU.

plotNotBoxPlot(data,num_cats,num_subcats,x_labels)



