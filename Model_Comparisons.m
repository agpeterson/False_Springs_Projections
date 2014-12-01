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
PATH_PREFIX = '/media/alexander/Vault/Bioclimate/';
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

% Create constant for number of years, lat, lon, and models.
N_LAT = 585;
N_LON = 1386;
N_MDL = 20;
N_VAR = 2;

% Create and preallocate writable matfile.
results = matfile('Results.mat','Writable',true);
results.gu_historic = NaN(N_LAT,N_LON,N_MDL);
results.lsf_historic = NaN(N_LAT,N_LON,N_MDL);

% Access MACAv2-METDATA lat/lon points, change to negative lon, and write to
% results file.
load('MACAv2METDATA_Coordinates.mat')
lon = lon - 360;
results.lat = lat;
results.lon = lon;

% Create years variables.
YRS_HST_IND = 1:56;
YRS_R45_IND = 57:150;
YRS_R85_IND = 151:244;


%%=============================================================================
% Main loop to access and extract model results.
%==============================================================================

% Iterate over variables and models.
for i=1:length(VARIABLE)
	for j=1:length(MODEL)

		% Concatenate to create filename for pointer; print 1 if exists.
		filename = [PATH_PREFIX,char(VARIABLE{i}),char(MODEL{j}),FILE_EXT]
		exist filename
	
		% Create pointer file.
		file = matfile(filename);

		% Switch on i to extract 1980-2009 values for each model and average 
		% for historic climatology.
		if i == 1			
			results.gu_historic(:,:,j) = squeeze(nanmean(file.gsi_CONUS...
														(YRS_HST_IND,:,:),1));
		else
			results.lsf_historic(:,:,j) = squeeze(nanmean(file.lsf_CONUS...
														(YRS_HST_IND,:,:),1));
		end

	end 	% j; 1:length(MODEL)

end 	% i; 1:length(VARIABLE)


%%=============================================================================
% Map climatologies.
%==============================================================================

% Initialize variables.
prj = 'Albers Equal-Area Conic';
min_val = 1;
max_val = 181;
val_step = 20;
lat_buffer = 2;
lon_buffer = 2;
cb_type = 'div';
cb_color = 'RdBu';
cb_units = 'Day of Year';
cb_flip = 'No Flip';

% Iterate over models for mapping.
figure('Position',[100 100 1000 618]);
for i=1:20
	
	% Extract data to plot and generate map title.
	data = results.gu_historic(:,:,i);
	map_title = [char(MODEL{i}),' Mean GU (1980-2009)'];

	% Call subaxis function and mapGriddedData function.
	figure('Position',[100 100 1000 618]);
	mapGriddedData(data,prj,min_val,max_val,val_step,...
				   lat,lon,lat_buffer,lon_buffer,...
				   map_title,cb_type,cb_color,cb_units,cb_flip)

end 	% i: 1:20





