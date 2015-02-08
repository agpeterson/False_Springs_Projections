%%===========================================================================
%
% NAME:   calcFSEIChangeSignif.m
% AUTHOR: Alexander Peterson
% DATE:   23 Jan. 2015
% DESC:   Determines changes and significance in FSEI between two sets of
%				false springs.
%
%============================================================================

function [fsei_change,fsei_significance] = calcFSEIChangeSignif(fs1,fs2)


%%============================================================================
% Initial variables and constants.
%-----------------------------------------------------------------------------

% Set parallel pool options.
opts = statset('UseParallel','always');

% Open parallel processor pool.
if matlabpool('size') == 0
    matlabpool open local 12
end

% Create N_* constants. 
N_YRS = size(fs1,1);
N_LAT = size(fs1,2);
N_LON = size(fs1,3);



%%============================================================================
% Bootstrap resampling.
%-----------------------------------------------------------------------------

% Iterate over lon and call bootstrp() function using 1000 samples, taking 
% the mean of the difference between the future and historic periods.
disp('Preallocating bootstrap variable...')
fs_bootstrap = NaN(1000,N_LAT,N_LON);

disp('Bootstrap function call...')
for lon=1:N_LON
    fs_bootstrap(:,:,lon) = bootstrp(1000,@nanmean,...
    	fs1(:,:,lon)-fs2(:,:,lon),'Options',opts);
end



%%============================================================================
% Change and significance.
%-----------------------------------------------------------------------------

% Magnitude of change.
disp('Calculating magnitude of change...')
fsei_change = squeeze(nanmean(fs_bootstrap,1));

% Find lat/lon greater/less than 0 to determine increasing/decreasing
% false springs.
disp('Determining statistical significance...')
fs_greater = squeeze(nansum(fs_bootstrap > 0));
fs_less = squeeze(nansum(fs_bootstrap < 0));
fs_increase = fs_greater >= 0.95*1000;
fs_decrease = fs_less >= 0.95*1000;

% Iterate over lat/lon and classify significance.
for x=1:N_LON
    for y=1:N_LAT
        if fs_decrease(y,x) == 1 || fs_increase(y,x) == 1
            fsei_significance(y,x) = 1;
        else
            fsei_significance(y,x) = 0;
        end
    end
end



end 	% Function.