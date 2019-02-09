function small_multiples_spectogram(signal,t,f,varargin)
%% DJC - 5-15-2018 small multiples spectrogram
% plot small mutliple plots - time x freq x channels x trials is input
% then it's averaged

signal = mean(signal,4);


% defaults
type1 = [];
type2 = [];
average = 0;
xlimits = [-200 1000];

for i=1:2:(length(varargin)-1)
    if ~ischar (varargin{i}),
        error (['Unknown type of optional parameter name (parameter' ...
            ' names must be strings).']);
    end
    % change the value of parameter
    switch lower (varargin{i})
        case 'type1'
            type1 = varargin{i+1};
        case 'type2'
            type2 = varargin{i+1};
        case 'average'
            average = varargin{i+1};
        case 'xlims'
            xlimits = varargin{i+1};
            
    end
end

%
totalFig = figure;
totalFig.Units = 'inches';
totalFig.Position = [   10.4097    3.4722   13.2708   10.4514];


p = numSubplots(size(signal,3));
%min_c = squeeze(min(min(min(signal))));
%max_c = squeeze(max(max(max(signal))));
minC = -3;
maxC = 3;
%cmap=flipud(cbrewer('div', 'RdBu', 13));
load('america');
cmap = cm;
CT = cm;
colormap(cmap)

for idx=1:size(signal,3)
    %smplot(p(1),p(2),idx,'axis','on')
    smplot(p(1),p(2),idx,'axis','on')
    
    if ismember(idx,type1)
        surf(1e3*t,f,zeros(size(signal(:,:,idx))),'edgecolor','k');
        title([num2str(idx)],'Color',CT(3,:))
        
    elseif ismember(idx,type2)
        surf(1e3*t,f,signal(:,:,idx),'edgecolor','k');
        title([num2str(idx)],'Color',CT(2,:))
    else
        surf(1e3*t,f,signal(:,:,idx),'edgecolor','k');
        
        title([num2str(idx)],'color','k')
    end
    view(0,90);
    axis tight;
    
    colormap(cmap);
    set_colormap_threshold(gcf, [-0.5 0.5], [minC maxC], [1 1 1])
    
    axis off
    axis tight
    hold on
    plot3([0,0],[0 300],[1000,1000],'k','linewidth',2)
    ylim([1 300])
    xlim(xlimits)
    
    %  vline(0,'k')
    xlabel('time (ms)');
    ylabel('frequency (Hz)');
    
end

colorbar()

end