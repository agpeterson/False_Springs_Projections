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
PATH_PREFIX = 'False_Springs_Data/';	% Thunder
%PATH_PREFIX = '/media/alexander/Vault/Bioclimate/';	% Local
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


%%=============================================================================
% Main loop to access file and derive metrics.
%==============================================================================

% Preallocate GU, LSF, and FSEI matfiles.
disp('Preallocating FS file...')
fsei_file = matfile('FSEI.mat','Writable',true);
fsei_file.clm_mean = NaN(N_CLM,N_LAT,N_LON,N_MDL,'single');

% Iterate over models for both variables.
for i=1:N_MDL

	% Concatenate to create filename for pointer.
	gu_filename = [PATH_PREFIX,char(VARIABLE{1}),char(MODEL{i}),FILE_EXT];
	lsf_filename = [PATH_PREFIX,char(VARIABLE{2}),char(MODEL{i}),FILE_EXT];
	
	% Create file pointers.
	gu_data = matfile(gu_filename);
	lsf_data = matfile(lsf_filename);

	% Print file to screen.
	disp(['Accessed file ', gu_filename])
	disp(['Accessed file ', lsf_filename])

	% Load model data into memory.
	gu = double(gu_data.gsi_CONUS);
	lsf = double(lsf_data.lsf_CONUS);


	%%=========================================================================
	% Derive false springs.
	%==========================================================================

	% Preallocate FS variable.
	disp('Preallocating FS variable...')
	fs = NaN(N_YRS,N_LAT,N_LON);

	% Iterate over lon, lat, and years.
	disp('Classifying false springs...')
	for j=1:N_LON
		for k=1:N_LAT
			for l=1:N_YRS
				if isnan(lsf(l,k,j)) || isnan(gu(l,k,j))
					fs(l,k,j) = NaN;
				elseif lsf(l,k,j) >= (7 + gu(l,k,j))
					fs(l,k,j) = 1;
				else
					fs(l,k,j) = 0;
				end 	% If statement.
			end 	% l; 1:N_YRS
		end 	% k; 1:N_LAT
	end 	% j; 1:N_LON

	% Preallocate FSEI variables.
	disp('Preallocating FSEI variables...')
	fs_sum = NaN(N_CLM,N_LAT,N_LON);
	gu_sum = NaN(N_CLM,N_LAT,N_LON);
	fsei = NaN(N_CLM,N_LAT,N_LON);

	% Derive FSEI.
	disp('Deriving FSEI...')
	for m=1:N_CLM
		fs_sum(m,:,:) = nansum(fs(clm_ind{m},:,:),1);	% Sums ones.
		gu_sum(m,:,:) = sum(~isnan(gu(clm_ind{m},:,:)),1);	% Sums # of elements.
		fsei(m,:,:) = (fs_sum(m,:,:) ./ gu_sum(m,:,:)) * 100;
	end

	% Write output to file.
	disp('Writing output to file...')
	fsei_file.clm_mean(:,:,:,i) = single(fsei);

	% Clear variables.
	clear gu_filename lsf_filename gu_data lsf_data fs fs_sum lsf_sum fsei
	clear lsf gu

end 	% j; 1:N_MDL


%%=============================================================================
% Derive climatological means and differences.
%==============================================================================

disp('Preallocating GU file...')
gu_file = matfile('GU.mat','Writable',true);
gu_file.clm_mean = NaN(N_CLM,N_LAT,N_LON,N_MDL,'single');
gu_file.clm_diff = NaN(N_CLM,N_LAT,N_LON,N_MDL,'single');

disp('Preallocating LSF file...')
lsf_file = matfile('LSF.mat','Writable',true);
lsf_file.clm_mean = NaN(N_CLM,N_LAT,N_LON,N_MDL,'single');
lsf_file.clm_diff = NaN(N_CLM,N_LAT,N_LON,N_MDL,'single');

% Iterate over models for both variables.
for i=1:N_MDL

	% Concatenate to create filename for pointer.
	gu_filename = [PATH_PREFIX,char(VARIABLE{1}),char(MODEL{i}),FILE_EXT];
	lsf_filename = [PATH_PREFIX,char(VARIABLE{2}),char(MODEL{i}),FILE_EXT];
	
	% Create file pointers.
	gu_data = matfile(gu_filename);
	lsf_data = matfile(lsf_filename);

	% Print file to screen.
	disp(['Accessed file ', gu_filename])
	disp(['Accessed file ', lsf_filename])

	% Load model data into memory.
	gu = double(gu_data.gsi_CONUS);
	lsf = double(lsf_data.lsf_CONUS);


	%%=========================================================================
	% Calculate climatological means and differences.
	%==========================================================================

	% Preallocate GU, LSF, and FSEI variables.
	disp('Preallocating GU and LSF variables...')
	gu_clm_mean = NaN(N_CLM,N_LAT,N_LON);
	gu_clm_diff = NaN(N_CLM,N_LAT,N_LON);
	lsf_clm_mean = NaN(N_CLM,N_LAT,N_LON);
	lsf_clm_diff = NaN(N_CLM,N_LAT,N_LON);

	% Iterate over N_CLM to aggregate.
	disp('Calculating climatological means and differences...')
	for j=1:N_CLM

		% Calculate climatological mean.
		gu_clm_mean(j,:,:) = nanmean(gu(clm_ind{j},:,:),1);
		lsf_clm_mean(j,:,:) = nanmean(lsf(clm_ind{j},:,:),1);

		% Take difference between historical and future periods.
		gu_clm_diff(j,:,:) = gu_clm_mean(j,:,:) - gu_clm_mean(1,:,:);
		lsf_clm_diff(j,:,:) = lsf_clm_mean(j,:,:) - lsf_clm_mean(1,:,:);
		
	end

	% Write output to file.
	disp('Writing output to files...')
	gu_file.clm_mean(:,:,:,i) = single(gu_clm_mean);
	gu_file.clm_diff(:,:,:,i) = single(gu_clm_diff);

	lsf_file.clm_mean(:,:,:,i) = single(lsf_clm_mean);
	lsf_file.clm_diff(:,:,:,i) = single(lsf_clm_diff);

	% Clear variables.
	clear gu_filename lsf_filename gu_data lsf_data lsf gu

end 	% j; 1:N_MDL


% Do the same as above but for FSEI differences.
disp('Preallocating FSEI file...')
fsei_file = matfile('FSEI.mat','Writable',true);
fsei_file.clm_diff = NaN(N_CLM,N_LAT,N_LON,N_MDL,'single');
fsei_file.clm_chng = NaN(N_CLM,N_LAT,N_LON,N_MDL,'single');

% Iterate over models.
for i=1:N_MDL

	% Load FSEI data into memory.
	disp('Loading model FSEI data into memory...')
	fsei = double(fsei_file.clm_mean(:,:,:,i));
	fsei_hst = squeeze(fsei(1,:,:));


	% Preallocate per model.
	disp('Preallocating variables per model...')
	fsei_diff = NaN(N_CLM,N_LAT,N_LON);
	fsei_chng = NaN(N_CLM,N_LAT,N_LON);

	% Iterate over N_CLM to aggregate.
	disp('Calculating FSEI differences and changes...')
	for j=1:N_CLM

		% Take difference between historical and future periods.
		fsei_diff(j,:,:) = squeeze(fsei(j,:,:)) - fsei_hst;

		% Calculate percent change between periods.
		fsei_chng(j,:,:) = ((squeeze(fsei(j,:,:)) - fsei_hst) ./ ...
			fsei_hst) * 100;

	end 	% j; 1:N_CLM

	% Write output to file.
	disp('Writing output to file...')
	fsei_file.clm_diff(:,:,:,i) = single(fsei_diff);
	fsei_file.clm_chng(:,:,:,i) = single(fsei_chng);

end 	% i; 1:N_MDL


%%=============================================================================
% Create historical climatology maps.
%==============================================================================

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

file = matfile('/media/alexander/Vault/Bioclimate/GU.mat')

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

% Export MACA coords to be used in Arc.
addpath(genpath('../Data'))
load MACAv2METDATA_Coordinates
lon = lon - 360;
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

% Regional historical means to test.
gu_hst = double(file.gu_mdl_mean(:,:,1));

for i=1:length(region_lat)
	data = gu_hst(region_lat{i},region_lon{i});
	data_rshp = reshape(data,[size(data,1)*size(data,2) 1]);
	test_mean(i,:) = nanmean(data_rshp,1);
end


