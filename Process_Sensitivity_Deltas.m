%%============================================================================
%
% NAME:     Process_Sensitivity_Deltas.m
% AUTH:     Alexander Peterson
% DATE:     8 Feb. 2015
% DESC:     This script accesses MACAv2-METDATA tmin files to create mean tmin
%           and standard deviation for use in sensitivity analysis of false 
%           springs.
%
%=============================================================================



%%============================================================================
% Initial variables and constants.
%-----------------------------------------------------------------------------

% Variable and model names.
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

% File prefix and suffix to be concatenated for file access.
FILE_PREFIX = '/storage/DOWNSCALED/CMIP5/MACAv2-METDATA/maca_v2_metdata_2var_10pat_CONUS_';
FILE_SUFFIX = {'_historical_tasmin.mat';...
               '_rcp85_tasmin.mat'};

% Break CONUS grid into regional subsets. The more powerful the machine, the
% larger the regional subsets can be.
LON_START = [1 701];
LON_END = [700 1386];

% Create iteration constants.
N_LAT = 585;
N_LON = 1386;
N_MDL = 20;
N_EXP = 2;


%-----------------------------------------------------------------------------
% Create and preallocate output file.
%-----------------------------------------------------------------------------
out_file = matfile('MACA_Sensitivity_Deltas.mat','Writable',true);
out_file.hst_mean = NaN(N_LAT,N_LON,N_MDL,'single');
out_file.hst_stnd = NaN(N_LAT,N_LON,N_MDL,'single');
out_file.fut_mean = NaN(N_LAT,N_LON,N_MDL,'single');
out_file.fut_stnd = NaN(N_LAT,N_LON,N_MDL,'single');



%%============================================================================
% Main loop to access matfiles. Splits on N_EXP to determine number of years
% (56 for historical; 30 for future).
%-----------------------------------------------------------------------------

disp('Entering model loop...')
for i=1:N_MDL

	disp('Entering experiment loop...')
	for j=1:N_EXP

		% Concatenate prefixes and suffixes to create file name.
		file_name = [FILE_PREFIX,char(MDL_NAME(i)),char(FILE_SUFFIX(j))];
		file = matfile(file_name);
		disp(['Accessed ',char(MDL_NAME(i)),char(FILE_SUFFIX(j))])


        %---------------------------------------------------------------------
        % Historic file access and processing. Create n_yrs variable and 
        % preallocate before entering LON loop.
        %---------------------------------------------------------------------
        if j == 1

            n_yrs = 56;

            disp('Preallocating stnd variable...')
            tmin_mean = NaN(n_yrs,N_LAT,N_LON,'single');
            tmin_stnd = NaN(n_yrs,N_LAT,N_LON,'single');
            tmin_mean_climo = NaN(N_LAT,N_LON,'single');
            tmin_stnd_climo = NaN(N_LAT,N_LON,'single');


            %-----------------------------------------------------------------
            % Iterate over LON to extract data. Pull days 60:151
            % (corresponding to Mar-May) for years 1-56 (1950-2005). Take the
            % seasonal mean and standard deviation, then calculate the
            % climatological mean of both.
            %-----------------------------------------------------------------
            disp('Entering LON loop...')
            for k=1:length(LON_START)

                % Break CONUS into regional lon subset.
                lon_subset = [LON_START(k):LON_END(k)];
                n_lon = length(lon_subset);

                % Extract data.
                disp('Loading historical TMIN data for days 60-151...')
                t_tmin = file.data(60:151,1:56,:,lon_subset) - 273.15;

                % Calculate mean and stnd.
                disp('Taking mean and stnd for all years...')
                t_mean = squeeze(nanmean(t_tmin,1));
                t_stnd = squeeze(nanstd(t_tmin));
                tmin_mean(:,:,lon_subset) = t_mean;
                tmin_stnd(:,:,lon_subset) = t_stnd;

                % Calculate climatological means.
                disp('Taking climatological mean...')
                t_mean_climo = squeeze(nanmean(t_mean,1));
                t_stnd_climo = squeeze(nanmean(t_stnd,1));
                tmin_mean_climo(:,lon_subset) = t_mean_climo;
                tmin_stnd_climo(:,lon_subset) = t_stnd_climo;

                % Clear variables.
                clear lon_subset n_lon t_tmin t_stnd t_climo

            end     % k; 1:length(LON_START)


            %-----------------------------------------------------------------
            % Write output to file.
            %-----------------------------------------------------------------
            disp('Writing output to file...')
            out_file.hst_mean(:,:,i) = tmin_mean_climo;
            out_file.hst_stnd(:,:,i) = tmin_stnd_climo;

        
        %---------------------------------------------------------------------
        % Future file access and processing; doesn't require LON loop. Similar
        % to LON loop above.
        %---------------------------------------------------------------------
        else
            
            % Create n_yrs variable based on experiment.
            n_yrs = 30;

            % Load tmin data for days 60-151 (March-May).
            disp('Loading future TMIN data for days 60-151...')
            t_tmin = file.data(60:151,35:64,:,:) - 273.15;

            % Calculate mean and stnd.
            disp('Taking mean and stnd for all years...')
            tmin_mean = squeeze(nanmean(t_tmin,1));
            tmin_stnd = squeeze(nanstd(t_tmin));

            % Calculate climatological means.
            disp('Taking climatological mean...')
            tmin_mean_climo = squeeze(nanmean(tmin_mean,1));
            tmin_stnd_climo = squeeze(nanmean(tmin_stnd,1));

            % Write output to file.
            disp('Writing output to file...')
            out_file.fut_mean(:,:,i) = tmin_mean_climo;
            out_file.fut_stnd(:,:,i) = tmin_stnd_climo;

        end     % if statement.

    end 	% j; 1:N_EXP
end 	% i; 1:N_MDL



%%============================================================================
% Create deltas.
%-----------------------------------------------------------------------------

% Load data into memory.
load('MACA_Sensitivity_Deltas.mat')

% Iterate over models to create deltas between historic and future values.
for i=1:N_MDL
    delta_mean(:,:,i) = fut_mean(:,:,i) - hst_mean(:,:,i);
    delta_stnd(:,:,i) = (fut_stnd(:,:,i) - hst_stnd(:,:,i)) ./ hst_stnd(:,:,i);
end

% Append to matfile.
save('MACA_Sensitivity_Deltas.mat','delta_mean','delta_stnd','-append')



%%============================================================================
% Exit program.
%-----------------------------------------------------------------------------
disp('Exiting MATLAB...')
exit