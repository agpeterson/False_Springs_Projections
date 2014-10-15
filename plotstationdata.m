%==============================================================================
% NAME:   plotstationdata.m
% DESC:   This function plots station data.
% IN:     Variable trend, significance, latitude, longitude, max value, and
%         min value.
% OUT:    Figure.
% CALL:   m_map.
% AUTH:   John Abatzoglou, heavily modified by Alexander Peterson. Last
%         modified 08/24/2014.
% NOTE:   Removed similar code from directory; needs further refactoring and
%         clean-up.
%==============================================================================

function plotstationdata(varTrend,varSig,lat,lon,maxVal,minVal);

% Set figure properties, projection.
set(gcf,'PaperPositionMode','auto');
m_proj('miller','long',[min(lon(:))-.1 max(lon(:))+.1],'lat',...
    [min(lat(:))-.1 max(lat(:))+.1]);
hold on;
load usahi

% Draw statelines.
for i=1:51
    m_plot(stateline(i).long,stateline(i).lat,'k')
end

% Colormap.
%colormap(flipud(jet(length(yTick)+1)));
colormap(flipud(jet(12)))
f=colormap;
j=f;

% Colorbar.
%v=colorbar('h');
%set(v,['y','TickLabel'],[yTick]);
%set(v,'ylim',([minVal maxVal]));

% Set bounds on data.
f=find(varTrend>maxVal);varTrend(f)=maxVal;
f=find(varTrend<minVal);varTrend(f)=minVal;

% Trend stuff (not sure).
diffVal=maxVal-minVal;
trendCeiling = (varTrend-minVal)/diffVal*11+1;
trendCeiling = ceil(trendCeiling);
%ranges=[1:1:length(yTick)+1];
%for i=1:length(varTrend)
%    for w=1:length(yTick)
%        if w==1
%            if varTrend(i)<yTick(w)
%                trendCeiling(i)=ranges(w);
%            end
%        elseif varTrend(i)>=(yTick(w-1)) && varTrend(i)<(yTick(w))
%            trendCeiling(i)=ranges(w);
%        elseif varTrend(i)>=(yTick(length(yTick)))
%            trendCeiling(i)=ranges(length(yTick)+1);
%        end
%    end
%end

% Plot stations.
for i=1:length(varTrend)
    if isnan(varTrend(i))==0
        lat2=lat(i);lon2=lon(i);
        if ~isnan(trendCeiling(i))
            h=m_plot(lon2,lat2,'s');
            set(h,'markersize',5,...
                'markerfacecolor',j(trendCeiling(i),:),...
                'markeredgecolor',j(trendCeiling(i),:));
            if varSig(i)==1
                set(h,'markersize',8,...
                    'markerfacecolor',j(trendCeiling(i),:),...
                    'markeredgecolor',j(trendCeiling(i),:),'linewidth',1,'marker','o');
            end
            hold on;
        end
    end
end

set(gca,'climmode','manual','clim',[minVal,maxVal]);

m_grid('linewi',0.0000001,'linest','none','color','k','tickdir','out',...
    'fontsize',16,'xtick',[],'ytick',[]);
set(findobj('tag','m_grid_color'),'facecolor','none')
set(gca,'fontsize',16);

% Colorbar.
v=colorbar('h');
set(v,'fontsize',16);
m_text(-99.5,18,'Difference (Days)','fontsize',18)

% Set position on screen.
set(gcf,'position',[ 50 50 1130 739]);

end