function small_multiples_fourier(P1,f,varargin)
% DJC - 2-18-2018 - small multiples plot for visualizing the results of
% time x channels x trials

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get inputs
p = inputParser;

addRequired(p,'P1',@isnumeric);
addRequired(p,'f',@isnumeric);

addParameter(p,'type1',[],@isnumeric);
addParameter(p,'type2',[],@isnumeric);
addParameter(p,'newFig',1,@(x) x==0 || x ==1)
addParameter(p,'plotLog',0,@(x) x==0 || x ==1)

p.parse(P1,f,varargin{:});

P1 = p.Results.P1;
f = p.Results.f;
type1 = p.Results.type1;
type2 = p.Results.type2;
newFig = p.Results.newFig;
plotLog = p.Results.plotLog;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define new figure
if newFig
    totalFig = figure;
    totalFig.Units = 'inches';
    totalFig.Position = [1 1 7 7];
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
p = vizFunc.numSubplots(size(P1,2));

nullSig = zeros(length(f),1);

if ~plotLog
    for idx=1:size(P1,2)
        %smplot(p(1),p(2),idx,'axis','on')
        plt_sub = vizFunc.smplot(p(1),p(2),idx,'axis','on');
        
        if ismember(idx,type1)
            plot(f,nullSig,'Color',CT(3,:),'LineWidth',2)
            title([num2str(idx)],'Color',CT(3,:))
        elseif ismember(idx,type2)
            plot(f,P1(:,idx),'Color',CT(2,:),'LineWidth',2)
            title([num2str(idx)],'Color',CT(2,:))
        else
            plot(f,P1(:,idx),'Color',CT(1,:),'LineWidth',2)
            title([num2str(idx)],'color',CT(1,:))
        end
        
        axis off
        axis tight
        xlim([0 500])
        ylim([0 2e-5])
        
    end
    
    %xlabel('f (Hz)')
    %ylabel('|P1(f)|')
    
    if ~newFig
        obj = vizFunc.scalebar;
        obj.XLen = 50;              %X-Length, 10.
        obj.XUnit = 'f (Hz)';            %X-Unit, 'm'.
        obj.YLen = 1e-5;
        obj.YUnit = 'Power';
        
        obj.Position = [5,0.25e-5];
        obj.hTextX_Pos = [5,0.5e-5]; %move only the LABEL position
        obj.hTextY_Pos =  [-10,0.5e-5];
        obj.hLineY(2).LineWidth = 5;
        obj.hLineY(1).LineWidth = 5;
        obj.hLineX(2).LineWidth = 5;
        obj.hLineX(1).LineWidth = 5;
        obj.Border = 'LL';          %'LL'(default), 'LR', 'UL', 'UR'
    end
else
    for idx=1:size(P1,2)
        %smplot(p(1),p(2),idx,'axis','on')
        plt_sub = vizFunc.smplot(p(1),p(2),idx,'axis','on');
        
        if ismember(idx,type1)
            semilogy(f,nullSig,'Color',CT(3,:),'LineWidth',2)
            title([num2str(idx)],'Color',CT(3,:))
        elseif ismember(idx,type2)
            semilogy(f,P1(:,idx),'Color',CT(2,:),'LineWidth',2)
            title([num2str(idx)],'Color',CT(2,:))
        else
            semilogy(f,P1(:,idx),'Color',CT(1,:),'LineWidth',2)
            title([num2str(idx)],'color',CT(1,:))
        end
        
        axis off
        axis tight
        xlim([0 1000])
        ylim([10e-14 10e-4])
        
    end
    
    %xlabel('f (Hz)')
    %ylabel('|P1(f)|')
    %
    %     obj = scalebar;
    %     obj.XLen = 50;              %X-Length, 10.
    %     obj.XUnit = 'f (Hz)';            %X-Unit, 'm'.
    %     obj.YLen = 1e-5;
    %     obj.YUnit = 'Power';
    %
    %     obj.Position = [5,0.25e-5];
    %     obj.hTextX_Pos = [5,0.5e-5]; %move only the LABEL position
    %     obj.hTextY_Pos =  [-10,0.5e-5];
    %     obj.hLineY(2).LineWidth = 5;
    %     obj.hLineY(1).LineWidth = 5;
    %     obj.hLineX(2).LineWidth = 5;
    %     obj.hLineX(1).LineWidth = 5;
    %     obj.Border = 'LL';          %'LL'(default), 'LR', 'UL', 'UR'
end



end