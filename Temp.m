%%=============================================================================
% Access LSF and GSI model output.
%==============================================================================

% Path suffix and prefix strings to be concatenated for model access.
FILE_PREFIX = 'False_Springs_Data/';	% Thunder
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

% Create path suffix and prefix strings to be concatenated.
TMIN_PREFIX = '/storage/DOWNSCALED/CMIP5/MACAv2-METDATA/maca_v2_metdata_2var_10pat_CONUS_';
TMIN_SUFFIX = {'_historical_tasmin';...
               '_rcp45_tasmin';...
               '_rcp85_tasmin'};

% Create constants for number of years, lat, lon, variables, and models.
N_LAT = 585;
N_LON = 1386;
N_MDL = 20;
N_VAR = 3;
N_YRS = 244;

% Break CONUS grid into regional subsets. The more powerful the machine, the
% larger the regional subsets can be.
LAT_START = [1 301];
LAT_END = [300 585];
LON_START = [1 501 1001];
LON_END = [500 1000 1386];

% Create index for historical and future years.
%HST_INDEX = 1:56;
%FUT_INDEX = [57:150] - 57 + 1;         % Sets FUT_INDEX as 1:94 for file access.

% Iterate over models.
for i=1:N_MDL

    % Create output file based on model name.
    out_filename = ['tmin_',char(MDL_NAME{i}),FILE_EXT];
    disp(['Creating file ', out_filename]);

    % Create matfile pointer and preallocate.
    out_file = matfile(out_filename,'Writable',true);
    out_file.tmin_avg_fut = NaN(30,N_LAT,N_LON,'single');
    out_file.tmin_min_fut = NaN(30,N_LAT,N_LON,'single');

	% Access GU/LSF file and load into memory.
	% Concatenate to create filename for pointer.
	gu_filename = [FILE_PREFIX,char(VAR_NAME{1}),char(MDL_NAME{i}),FILE_EXT];
    disp(['Accessing file ', gu_filename]);

	% Create file pointer and load into memory.
	gu_file = matfile(gu_filename);
	gu = gu_file.gsi_CONUS;

	% Find dates +/-30 days from GU/LSF for each year.
    disp('Finding +/-21 days from annual GU date...');
	doy = NaN(N_LAT,N_LON,'single');
	day_start = NaN(N_YRS,N_LAT,N_LON,'single');
	day_end = NaN(N_YRS,N_LAT,N_LON,'single');
	for j=1:N_YRS
		doy = squeeze(gu(j,:,:));
		day_start(j,:,:) = doy - 21;
		day_end(j,:,:) = doy + 21;
	end
    day_start(day_start < 1) = 1;
    day_end(day_end > 192) = 192;

	% Access tmin file and load above dates into memory.
	for k=1:N_VAR

		% Create path string for each file and set pointer to matfile.
        tmin_filename = [TMIN_PREFIX,char(MDL_NAME{i}),char(TMIN_SUFFIX{k}),FILE_EXT];
        tmin_file = matfile(tmin_filename);
        disp(['Accessing file ', tmin_filename]);

        % Create variable to hold number of years, switching on experiment.
        if k == 1
            yr_index = [1:56];
        elseif k == 2
            yr_index = [215:244];
        else k == 3
            yr_index = [215:244];
            % yr_index = [151:244];
     	end;
        n_yrs = length(yr_index);

        % Iterate over lons.
        disp('Iterating over regions...')
        for x=1:length(LON_START)
        
            % Break CONUS into regional lon subset.
            lon_subset = [LON_START(x):LON_END(x)];
            n_lon = length(lon_subset);
        
            for y=1:length(LAT_START)

                [x y]
            
            	% Break CONUS into regional lat, day, and vpd subsets.
                lat_subset = [LAT_START(y):LAT_END(y)];
                n_lat = length(lat_subset);

                % Create temporary variable to store daily and yearly data for
                % each lat/lon subset.
                disp('Preallocating temporary tmin variables...')
                t_tmin_avg = NaN(n_yrs,n_lat,n_lon,'single');
                t_tmin_min = NaN(n_yrs,n_lat,n_lon,'single');

                % Load tmin data into memory and convert from K to C.
                disp('Loading tmin data...')
                if k == 1
                    t_var = double(tmin_file.data(1:192,1:56,...
                        lat_subset,lon_subset)) - 273.15;
                else
                    t_var = double(tmin_file.data(1:192,66:95,...
                        lat_subset,lon_subset)) - 273.15;
                end

                % Get day range for lat/lon/year.
                disp('Calculating mean mean and and mean min tmin...')
                for l=1:n_lon
                	for m=1:n_lat
                		for n=1:n_yrs

                			% Create day range.
                			days = day_start(yr_index(n),lat_subset(m),lon_subset(l)):...
                                   day_end(yr_index(n),lat_subset(m),lon_subset(l));

                			% Find yearly avg/min.
                            if isnan(days)
                                t_tmin_avg(n,m,l) = NaN;
                                t_tmin_min(n,m,l) = NaN;
                            else
                                t_tmin_avg(n,m,l) = nanmean(t_var(days,n,m,l),1);
                                t_tmin_min(n,m,l) = min(t_var(days,n,m,l));
                            end     % if statement

                		end 	% n; 1:n_yrs
                	end 	% m; 1:n_lat
                end 	% l; 1:n_lon

                % Write regional subsets to file.
                disp('Writing output to file...')
                out_file.tmin_avg_fut(1:30,lat_subset,lon_subset) = t_tmin_avg;
                out_file.tmin_min_fut(1:30,lat_subset,lon_subset) = t_tmin_min;

            end 	% y; 1:length(LAT_START)

        end     % x; 1:length(LON_START)

    end     % k; 1:N_VAR
end     % i; 1:N_MDL


% Average min tmins over 30 years; average tmin average 30 years.

figure()
scatter(reshape(min_diff,[585*1386 1]),reshape(avg_diff,[585*1386 1]));

figure()
scatter(reshape(avg_diff,[300*500 1]),reshape(lsf,[300*500 1]));


lsf_file = matfile('LSF.mat')
lsf = squeeze(lsf_file.clm_diff(5,1:300,1:500,1));













