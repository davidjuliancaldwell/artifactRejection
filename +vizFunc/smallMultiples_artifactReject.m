function smallMultiples_artifactReject(signal,t,varargin)
% DJC - 2-18-2018 - small multiples plot for visualizing the results of 
% time x channels x trials 

% defaults
type1 = [];
type2 = [];
average = 0;
newfig = true;
highlight_range = [];
xlims = [-200 1000];
ylims = [-200 200];

% variable input 
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
        case 'newfig'
            newfig = varargin{i+1};
        case 'highlight_range'
            highlight_range = varargin{i+1};
        case 'xlims'
            xlims = varargin{i+1};
        case 'ylims'
            ylims = varargin{i+1};
            
    end
end

%% define new figure 
if newfig
    totalFig = figure;
    totalFig.Units = 'inches';
    totalFig.Position = [   10.4097    3.4722   13.2708   10.4514];
    CT = vizFunc.cbrewer('qual','Accent',8);
    CT = flipud(CT);
else
    gcf;
    hold on
    CT = vizFunc.cbrewer('qual','Accent',8);
    CT = flipud(CT);
    CT(1,:) = CT(4,:);
end


% determine number of subplots using subplots helper function 
p = vizFunc.numSubplots(size(signal,2));

nullSig = zeros(length(t),1);

for idx=1:size(signal,2)
    %smplot(p(1),p(2),idx,'axis','on')
    plt_sub = vizFunc.smplot(8,8,idx,'axis','on');
    
    if average
        if ismember(idx,type1)
            plot(1e3*t,nullSig,'Color',CT(3,:),'LineWidth',2)
            title([num2str(idx)],'Color',CT(3,:))
        elseif ismember(idx,type2)
            plot(1e3*t,1e6*signal(:,idx),'Color',CT(2,:),'LineWidth',2)
            title([num2str(idx)],'Color',CT(2,:))
        else
            plot(1e3*t,1e6*signal(:,idx),'Color',CT(1,:),'LineWidth',2)
            title([num2str(idx)],'color',CT(1,:))
        end
        
    elseif ~average
        if ismember(idx,type1)
            plot(1e3*t,nullSig,'Color',CT(3,:),'LineWidth',2)
            title([num2str(idx)],'Color',CT(3,:))
        elseif ismember(idx,type2)
            plot(1e3*t,1e6*squeeze(signal(:,idx)),'Color',CT(2,:),'LineWidth',2)
            title([num2str(idx)],'Color',CT(2,:))
        else
            plot(1e3*t,1e6*squeeze(signal(:,idx)),'Color',CT(1,:),'LineWidth',2)
            title([num2str(idx)],'color',CT(1,:))
        end
        
    end
    
    axis off
    axis tight
    %xlim([-10 200])
    xlim(xlims)
    
    ylim(ylims)
    vizFunc.vline(0)
    
    if ~isempty(highlight_range)
        hColor = [116/255 255/255 112/255];
        y_range = ([-150 150]);
         vizFunc.highlight(plt_sub, highlight_range, y_range, hColor);
    end
       
    
end

obj = vizFunc.scalebar;
obj.XLen = 100;              %X-Length, 10.
obj.XUnit = 'ms';            %X-Unit, 'm'.
obj.YLen = 200;
obj.YUnit = '\muV';

obj.Position = [20,-130];
obj.hTextX_Pos = [5,-50]; %move only the LABEL position
obj.hTextY_Pos =  [-45,-40];
obj.hLineY(2).LineWidth = 5;
obj.hLineY(1).LineWidth = 5;
obj.hLineX(2).LineWidth = 5;
obj.hLineX(1).LineWidth = 5;
obj.Border = 'LL';          %'LL'(default), 'LR', 'UL', 'UR'

end