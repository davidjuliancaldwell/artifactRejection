function small_multiples_spectrogram(signal,t,f,varargin)
%% DJC - 5-15-2018 small multiples spectrogram
% plot small mutliple plots - time x freq x channels x trials is input
% then it's averaged

signal = mean(signal,4);

% defaults
type1 = [];
type2 = [];
average = 0;
xlims = [-200 1000];
ylims = [1 300];

for i=1:2:(length(varargin)-1)
    if ~ischar (varargin{i})
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
            xlims = varargin{i+1};
        case 'ylims'
            ylims = varargin{i+1};
            
    end
end

%
totalFig = figure;
totalFig.Units = 'inches';
totalFig.Position = [ 1.7708 1.5208 12.7917 8.2604];

p = numSubplots(size(signal,3));
%min_c = squeeze(min(min(min(signal))));
%max_c = squeeze(max(max(max(signal))));
minC = -3;
maxC = 3;
%cmap=flipud(cbrewer('div', 'RdBu', 13));
load('america');
cmap = cm;
CT = cm;

for idx=1:size(signal,3)
    %smplot(p(1),p(2),idx,'axis','on')
    smplot(p(1),p(2),idx,'axis','on')
    
    if ismember(idx,type1)
        s = surf(1e3*t,f,zeros(size(signal(:,:,idx))),'edgecolor','none');
        title([num2str(idx)],'Color',CT(3,:))
        
    elseif ismember(idx,type2)
        s = surf(1e3*t,f,signal(:,:,idx),'edgecolor','none');
        title([num2str(idx)],'Color',CT(2,:))
    else
        s =  surf(1e3*t,f,signal(:,:,idx),'edgecolor','none');
        title([num2str(idx)],'color','k')
    end
    
    %     %%Create vectors out of surface's XData and YData
    %     x=s.XData;
    %     y=s.YData;
    %     z=s.ZData;
    %     x=x(1,:);
    %     y=y(1,:);
    %     interestX = [xlims(1) xlims(end)];
    %     interestY = [ylims(1) ylims(end)];
    %     % Divide the lengths by the number of lines needed
    %
    %     % Plot the mesh lines
    %     % Plotting lines in the X-Z plane
    %     hold on
    %     for i = 1:2
    %         Y1 = interestY(i)*ones(size(x)); % a constant vector
    %         Z1 = zeros(size(x));
    %         plot3(x,Y1,Z1,'-k');
    %     end
    %     % Plotting lines in the Y-Z plane
    %     for i = 1:2
    %         X2 = interestX(i)*ones(size(y)); % a constant vector
    %         Z1 = zeros(size(X2));
    %         plot3(X2,y,Z1,'-k');
    %     end
    
    view(0,90);
    axis tight;
    xlim(xlims);
    ylim(ylims);
    colormap(cmap);
    vizFunc.set_colormap_threshold(gcf, [-0.5 0.5], [minC maxC], [1 1 1])
    axis off
    vline(0,'k')

end

colorbar()
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
end