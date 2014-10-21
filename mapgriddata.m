%%=============================================================================
% NAME:   mapgriddata.m
% AUTHOR: Alexander Peterson
% DATE:   30 Sept. 2014
% DESCR:  This script contains prototype code to plot the various components of
%		  false springs.
% IN:     
% OUT:    
% CALLS:  m_map; conus_grid.mat; usahi.mat; m_proj; m_coast; m_pcolor; m_grid;
%		  clbmap
%==============================================================================

function mapgriddata(data,min_val,max_val,val_step,lat,lon,...
						  map_title,cb_type,cb_color,cb_units,cb_flip)

% Add necessary paths and load data for plotting.
load('usahi.mat')

% Set variables for colorbar ticks and classes.
cb_ticks = [min_val:val_step:max_val];
n_classes = length(cb_ticks)-1;

% Create figure and set background to white.
fig = figure('Position',[100 100 1000 618]);	% Golden ratio.
set(gcf,'Color','W');

% Set projection and bounds.
m_proj('Albers Equal-Area Conic',...
       'lon',[min(lon(:))-2 max(lon(:))+2],...
       'lat',[min(lat(:))-2 max(lat(:))+2]);
hold on;

% Load and draw coastline.
m_coast('Linewidth',.1,'Color','k');

% Plot data.
m_pcolor(lon,lat,data);
shading flat;

% Remove default grid.
m_grid('xTick',[],'yTick',[]);

% Set colormap.
set(gca,'Climmode','Manual','Clim',[min_val max_val]);
if strcmp(cb_flip,'Flip') == 1
	colormap(flipud(cbrewer(cb_type,cb_color,n_classes)));
else
	colormap(cbrewer(cb_type,cb_color,n_classes));
end

% Draw statelines.
for i=1:51
    m_plot(stateline(i).long,stateline(i).lat,'k');
end

% Colorbar.
v = colorbar('h');

% Set colorbar label.
set(get(v,'yLabel'),'String',cb_units,'FontSize',16);

% Set color axes.
caxis([min(cb_ticks) max(cb_ticks)]);
set(v,'xTick',cb_ticks,'FontSize',14);

% Set title
title({''; map_title},'FontSize',18)












end