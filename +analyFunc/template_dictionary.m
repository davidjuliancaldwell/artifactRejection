function [processedSig,templateArrayCellOutput] = template_dictionary(templateArrayCell,templateTrial,rawSig,fs,varargin)
% Usage:  [processedSig,templateArrayCellOutput] = template_dictionary(templateArrayCell,templateTrial,rawSig,varargin)
%
% This function implements the template dictionary method. Briefly, the
% beginning and ending indices of stimulation artifact periods are
% extracted on a channel and trial-wise basis. From here,
%
% Arguments:
%   Required:
%   templateArrayCell - a cell array with a collection of templates on a
%   channel basis. Each entry, (e.g. templateArrayCell{channel}) has all of
%   the channelwise artifact epochs.
%   templateTrial - a cell array with a collection of templates on a
%   channel and trial basis. Each entry, (e.g.
%   templateTrial{channel}{trial} has all individual artifact epochs for a
%   given channel and trial.
%  rawSig - time x channels x trials raw signal
%          fs - sampling rate (Hz)
%
%
%   Optional:
%useFixedEnd - use a fixed end distance (1), or dynamically calculate the
%              offset of each stimulus pulse
%      fixedDistance - the maximum distance in ms to either look beyond
%        pre - the number of ms before which the stimulation pulse onset as
%              detected by a thresholding method should still be considered
%              as artifact
%       post - the number of ms after which the stimulation pulse onset as
%              detected by a thresholding method should still be
%              considered as artifact
%      plotIt - plot intermediate steps if true
%     goodVec - vector (e.g. [1 2 4 10]) of the channel indices to use. If
%              there are bad channels or stimulation channels for instance,
%              they should be excluded from goodVec
%  recoverExp - 1 or 0. If "1", try to recover an early cortical response
%  on a decaying exponential within the artifact transient.
% distanceMetricDbscan - distance metric to use with the DBScan dictionary
% building method.
%
% Returns:
%      startInds - cell array of the start indices each artifact for each
%      channel and trial - startInds{trial}{channel}
%       endsInds - cell array of the end indices of each artifact for each
%      channel and
%
%
% Copyright (c) 2018 Updated by David Caldwell
% University of Washington
% djcald at uw . edu
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;

addRequired(p,'templateArrayCell',@iscell);
addRequired(p,'templateTrial',@iscell);
addRequired(p,'rawSig',@isnumeric);
addRequired(p,'fs',@isnumeric);

addParameter(p,'plotIt',0,@(x) x==0 || x ==1);

addParameter(p,'distanceMetricDbscan','eucl',@isstr);
addParameter(p,'distanceMetricSigMatch','eucl',@isstr);

addParameter(p,'goodVec',[1:64],@isnumeric);
addParameter(p,'startInds',[],@iscell);
addParameter(p,'endInds',[],@iscell);

addParameter(p,'recoverExp',1,@(x) x==0 || x ==1);
addParameter(p,'maxAmps',ones(size(rawSig,2),size(rawSig,3)),@isnumeric)
addParameter(p,'maxLocation',15,@isnumeric);

addParameter(p,'normalize','preAverage',@isstr);
addParameter(p,'amntPreAverage',3,@isnumeric);
addParameter(p,'bracketRange',[-6:6],@isnumeric);

addParameter(p,'expThreshVoltageCut',95,@isnumeric);
addParameter(p,'expThreshDiffCut',95,@isnumeric);
addParameter(p,'chanInt',1,@isnumeric);

addParameter(p,'minPts',2,@isnumeric);
addParameter(p,'minClustSize',3,@isnumeric);
addParameter(p,'outlierThresh',0.95,@isnumeric);

p.parse(templateArrayCell,templateTrial,rawSig,fs,varargin{:});

templateArrayCell = p.Results.templateArrayCell;
templateTrial = p.Results.templateTrial;
rawSig = p.Results.rawSig;
fs =  p.Results.fs;

plotIt = p.Results.plotIt;
distanceMetricDbscan = p.Results.distanceMetricDbscan;
distanceMetricSigMatch = p.Results.distanceMetricSigMatch;
goodVec = p.Results.goodVec;
startInds = p.Results.startInds;
endInds = p.Results.endInds;

recoverExp = p.Results.recoverExp;
maxAmps = p.Results.maxAmps;

normalize = p.Results.normalize;
amntPreAverage = p.Results.amntPreAverage;
maxLocation = p.Results.maxLocation;
bracketRange = p.Results.bracketRange;

expThreshVoltageCut = p.Results.expThreshVoltageCut;
expThreshDiffCut = p.Results.expThreshDiffCut;

chanInt = p.Results.chanInt;

minPts = p.Results.minPts;
minClustSize = p.Results.minClustSize;
outlierThresh = p.Results.outlierThresh;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
templateArrayCellOutput = {};
processedSig = zeros(size(rawSig));
templateSubtractCell = {};
templateListVec = {};

fprintf(['-------Dictionary-------- \n'])

plotIt = 1;

for chan = goodVec
    fprintf(['-------Artifact Channel ' num2str(chan) ' -------- \n'])
    
    templateArray = templateArrayCell{chan};
    
    % extract max amplitude for a given channel
    maxAmpsChan = max(maxAmps(chan,:));
    
    % shorten data to be centered around the peak +/- the bracketRange. In
    % this way there is less clustering around non-discriminative data
    % points.
    templateArrayShortened = templateArray(maxLocation+bracketRange,:);
    [~,maxSub] = max(templateArrayShortened,[],1);
    maxSub = round(median(maxSub));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % data is in "templateArrayShorted". We will initiate a new HDBSCAN instance
    clusterer = HDBSCAN.HDBSCAN( templateArrayShortened');
    
    % we can view our data matrix size
    fprintf( 'Number of points: %i\n',clusterer.nPoints );
    fprintf( 'Number of dimensions: %i\n',clusterer.nDims );
    
    try
        % (1) directly set the parameters
        %         clusterer.minpts = 2;
        %         clusterer.minclustsize = 3;
        %         clusterer.outlierThresh = 0.95;
        clusterer.minpts = minPts;
        clusterer.minclustsize = minClustSize;
        clusterer.outlierThresh = outlierThresh;
        
        
        clusterer.metric = distanceMetricDbscan;
        clusterer.fit_model(); 			% trains a cluster hierarchy
    catch
        % (1) directly set the parameters
        clusterer.minpts = minPts+1;
        clusterer.minclustsize = minClustSize+1;
        clusterer.outlierThresh = outlierThresh;
        clusterer.metric = distanceMetricDbscan;
        clusterer.fit_model(); 			% trains a cluster hierarchy
    end
    clusterer.get_best_clusters(); 	% finds the optimal "flat" clustering scheme
    clusterer.get_membership();		% assigns cluster labels to the points in X
    
    % (2) call run_hdbscan() with optional inputs. This is the prefered/easier method
    %  clusterer.run_hdbscan( 10,20,[],0.85 );
    
    % Let's visualize the condensed cluster tree (the tree without spurious clusters)
    
    labels = clusterer.labels;
    vectorUniq = unique(labels);
    templateArrayExtracted = [];
    for i = vectorUniq'
        if i~=0
            meanTempArray = mean(templateArray(:,labels==i),2);
            templateArrayExtracted = [templateArrayExtracted (meanTempArray )]; %no subtraction
        end
    end
    
    %   if plotIt
    if plotIt && chan == chanInt
        %%
        colorFirst = [117 112 179;
            231 41 138;
            213 255 0;
            0 244 255;
            0.7 0.7 0.7]./255;
        
        figure
        currentFig = gcf;
        currentFig.Units = "inches";
        currentFig.Position = [1 1 1.8 1.8];
        set(gca,'fontsize',12);
        h = scatter(1e3*clusterer.data(:,maxSub),1e3*clusterer.data(:,maxSub+1),20,'filled');
        ylabel('V(t=1)')
        xlabel('V(t=2)')
        set( h.Parent,'tickdir','out','box','off' );
        
        tempLabels = clusterer.labels;
        tempLabels(tempLabels == 0) = max(unique(tempLabels))+1;
        [C,ia,ic] = unique(tempLabels);
        labelInfo.C{chan} = C;
        labelInfo.ic{chan} = ia;
        labelInfo.tempLabels{chan} = tempLabels;
        h.CData = ic;
        
        %minLength1 = 8;
        minLength1 = 5;
        minLength2 = length(unique(tempLabels))+2;
        CT = [colorFirst; cbrewerHelper.cbrewer('qual', 'Set1',max(minLength1,minLength2))];
        % colormapOrder = [3 4 6 8];
        % lengthCTvec = 1:size(CT,1);
        %   isMemberCT = ismember(lengthCTvec,colormapOrder);
        %   colormapOrder = [colormapOrder lengthCTvec(~isMemberCT)];
        %  CT = CT(colormapOrder,:);
        colors=CT;
        labelInfo.colors{chan} = colors;
        colorsNew = colors(1:length(unique(ic)),:);
        colormap(colorsNew);
        
        totalFig = figure;
        totalFig.Units = 'inches';
        totalFig.Position = [1 1 6.0417 9.4792];
        subplot(2,1,1)
        h = scatter(1e3*clusterer.data(:,maxSub),1e3*clusterer.data(:,maxSub+1),20,'filled');
        ylabel('time point 1 : voltage (mV)')
        xlabel('time point 2 : voltage (mV)')
        set( h.Parent,'tickdir','out','box','off' );
        h.CData = ic;
        colorsNew = colors(1:length(unique(ic)),:);
        colormap(colorsNew);
        
        set(gca,'fontsize',18)
        
        subplot(2,1,2)
        t = 1e3*[0:size(templateArrayShortened,1)-1]/fs;
        plot(t,1e3*templateArrayShortened,'color',[0.75 0.75 0.75])
        hold on
        h2 = plot(t,1e3*templateArrayExtracted(maxLocation+bracketRange,:),'linewidth',2);
        set(h2, {'color'}, num2cell(colors(1:length(C)-1,:),2));
        
        colormap( colors );
        ylabel('Voltage (mV)')
        xlabel('Time (ms)')
        set(gca,'fontsize',18)
        title('Dictionary of Templates from Raw Traces')
        
        index = 1;
        for labelsIndiv = 1:length(C)-1
            figure
            t = 1e3*[0:size(templateArrayShortened,1)-1]/fs;
            plot(t,1e3*templateArrayShortened(:,ic==labelsIndiv),'color',[0.75 0.75 0.75])
            hold on
            h2 = plot(t,1e3*templateArrayExtracted(maxLocation+bracketRange,index),'linewidth',2);
            set(h2, {'color'}, num2cell(colors(labelsIndiv,:),2));
            ylabel('Voltage (mV)')
            xlabel('Time (ms)')
            set(gca,'fontsize',18)
            title('Dictionary of Templates from Raw Traces')
            index = index + 1;
        end
        
        figure
        totalFig.Units = 'inches';
        totalFig.Position = [1 1 6.0417 9.4792];
        t = 1e3*[0:size(templateArrayShortened,1)-1]/fs;
        plot(t,1e3*templateArrayShortened,'color',[0.75 0.75 0.75])
        hold on
        h2 = plot(t,1e3*templateArrayExtracted(maxLocation+bracketRange,:),'linewidth',2);
        set(h2, {'color'}, num2cell(colors(1:length(C)-1,:),2));
        
        colormap( colors );
        
        ylabel('Voltage (mV)')
        xlabel('Time (ms)')
        set(gca,'fontsize',18)
        title('Dictionary of Templates from Raw Traces')
        clusterInt = clusterer;
    end
    templateListVec{chan} = templateArrayShortened;
    
    % assign templates to channel
    templateArrayCellOutput{chan} = templateArrayExtracted;
    
end

fprintf(['-------Finished clustering artifacts-------- \n'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% now do the template subtraction
for trial = 1:size(rawSig,3)
    
    if trial > 1
        firstLoopTrial = 0;
    else
        firstLoopTrial = 1;
    end
    
    for chan = goodVec
        
        firstLoopChan = 1;
        rawSigTemp = rawSig(:,chan,trial);
        templates = templateArrayCellOutput{chan};
        
        % add on the trial one
        %templates = [templates mean(templateTrial{chan}{trial},2)];
        templates = templates;
        % ensure no subtraction of exponential
        if recoverExp
            templates = analyFunc.recover_EP(templates,fs,'threshDiffCut',expThreshDiffCut,'threshVoltageCut',expThreshVoltageCut);
        end
        
        for sts = 1:length(startInds{trial}{chan})
            win = startInds{trial}{chan}(sts):endInds{trial}{chan}(sts);
            extractedSig = rawSigTemp(win);
            
            switch normalize
                case 'preAverage'
                    extractedSig = extractedSig - mean(extractedSig(1:amntPreAverage));
                case 'none'
                    extractedSig = extractedSig ;
                case 'firstSamp'
                    extractedSig = extractedSig - extractedSig(1);
                case 'mean'
                    extractedSig = extractedSig - mean(extractedSig);
            end
            
            % find best artifact
            templatesSts = templates(1:length(extractedSig),:);
            
            % make sure the bracket range does not exceed the first or last
            % sample of the template array
            
            sizeTemplates = size(templatesSts);
            bracketRange =  p.Results.bracketRange;
            bracketRangeMin = maxLocation+bracketRange(1);
            bracketRangeMax = maxLocation+bracketRange(end);
            
            adjustTemplates = false;
            if (bracketRangeMin < 1)
                bracketRangeMin = 0;
                adjustTemplates = true;
            end
            
            if (bracketRangeMax >= sizeTemplates(1))
                bracketRangeMax = sizeTemplates(1) - maxLocation;
                adjustTemplates = true;
            end
            
            if adjustTemplates
                bracketRange = [bracketRangeMin:bracketRangeMax];
            end
            
            if isempty(bracketRange)
               bracketRange = 1:length(templatesSts)-1;
               maxLocation = 1;
            end
            
            templatesStsShortened = templatesSts(maxLocation+bracketRange,:);
            extractedSigShortened = extractedSig(maxLocation+bracketRange,:);
            
            switch distanceMetricSigMatch
                case 'corr'
                    % correlation
                    correlations = corr(extractedSigShortened,templatesStsShortened);
                    [~,index] = max((correlations));
                    
                case 'cosine'
                    % cosine similarity
                    denominator = sqrt(sum(extractedSigShortened.*extractedSigShortened)).*sqrt(sum(templatesStsShortened.*templatesStsShortened));
                    numerator = (extractedSigShortened'*templatesStsShortened);
                    correlations = numerator./denominator;
                    [~,index] = max(abs(correlations));
                    
                case 'eucl'
                    % distance
                    v = templatesStsShortened - repmat(extractedSigShortened,1,size(templatesStsShortened,2));
                    distance = sum(v.*v);
                    [~,index] = min(distance);
                case 'dtw'
                    % dynamic time warping
                    sizeTemplates = size(templatesStsShortened,2);
                    dtwMat = zeros(1,sizeTemplates);
                    for index = 1:sizeTemplates
                        dtwMat(index) = dtw(templatesStsShortened(:,index),extractedSigShortened);
                    end
                    [~,index] = min(dtwMat);
            end
            
            templateSubtract = templatesSts(:,index);
            % which template best matched
            if firstLoopChan && firstLoopTrial
                templateSubtractCell{chan} = index;
            else
                templateSubtractCell{chan} = [templateSubtractCell{chan}; index];
            end
            
            
            % scaling
            scaling = mean([max(rawSigTemp(win))/max(templateSubtract),min(rawSigTemp(win))/min(templateSubtract)]);
            %    scaling = max(rawSigTemp(win))/max(templateSubtract);
            templateSubtract = templateSubtract*scaling;
            
            rawSigTemp(win) = rawSigTemp(win) - templateSubtract;
            
            if plotIt && chan == chanInt && (sts == 1 || sts == 2 || sts == 10) && firstLoopChan
                figure
                currentFig = gcf;
                currentFig.Units = "inches";
                currentFig.Position = [1 1 2 2];
                set(gca,'fontsize',8);
                t = [0:length(extractedSig)-1]/fs;
                
                plot(1e3*t,1e3*extractedSig,'linewidth',2,'color','k')
                hold on
                plot(1e3*t,1e3*templateSubtract,'linewidth',2,'color',[0 255 134]/255)
                plot(1e3*t,1e3*extractedSig-1e3*templateSubtract,'linewidth',2,'color',[204 85 0]/255)
                legend('Raw Signal','Template Selected','Recovered Signal');
                set(gca,'fontsize',8)
                title('Raw Signal, Template, and Recovered Signal')
                ylabel('Voltage (mV)')
                xlabel('Time (ms)')
            end
            
            firstLoopChan = 0;
            
            
        end
        
        if plotIt && (trial == 10 || trial == 15 || trial == 20) && chan == chanInt
            %%
            figure
            currentFig = gcf;
            currentFig.Units = "inches";
            currentFig.Position = [1 1 2.8330 3.6670];
            set(gca,'fontsize',12);
            hold on
            for indexPlot = 1:length(startInds{trial}{chan})
                tempBox = vizFunc.highlight(gca, [1e3*startInds{trial}{chan}(indexPlot)/fs 1e3*endInds{trial}{chan}(indexPlot)/fs], [-1e3*max(abs(rawSig(:,chan,trial))) 1e3*max(abs(rawSig(:,chan,trial)))], [.7 .7 .7]);
            end
            t = [0:length(rawSigTemp)-1]/fs;
            plot1 = plot(1e3*t,1e3*rawSig(:,chan,trial),'linewidth',2,'color','k');
            plot2 = plot(1e3*t',1e3*rawSigTemp,'linewidth',2,'color',[204 85 0]/255);
            xlim([900 1500])
            ylim([-1e3*max(abs(rawSig(:,chan,trial))) 1e3*max(abs(rawSig(:,chan,trial)))])
            legend([tempBox,plot1,plot2],{'Artifact Windows','Raw Signal','Recovered Signal'})
            title('Raw vs. Recovered Signal')
            xlabel('Time (ms)')
            ylabel('Voltage (mV)')
            
            figure
            currentFig = gcf;
            currentFig.Units = "inches";
            currentFig.Position = [1 1 2.8330 3.6670];
            set(gca,'fontsize',12);
            hold on
            for indexPlot = 1:length(startInds{trial}{chan})
                tempBox = vizFunc.highlight(gca, [1e3*startInds{trial}{chan}(indexPlot)/fs 1e3*endInds{trial}{chan}(indexPlot)/fs], [-1e3*max(abs(rawSig(:,chan,trial))) 1e3*max(abs(rawSig(:,chan,trial)))], [.7 .7 .7]);
            end
            t = [0:length(rawSigTemp)-1]/fs;
            plot1 = plot(1e3*t,1e3*rawSig(:,chan,trial),'linewidth',2,'color','k');
            plot2 = plot(1e3*t',1e3*rawSigTemp,'linewidth',2,'color',[204 85 0]/255);
            xlim([1000 1015])
            ylim([-1e3*max(abs(rawSig(:,chan,trial))) 1e3*max(abs(rawSig(:,chan,trial)))])
            legend([tempBox,plot1,plot2],{'Artifact Windows','Raw Signal','Recovered Signal'})
            title('Raw vs. Recovered Signal')
            xlabel('Time (ms)')
            ylabel('Voltage (mV)')
            
            figure
            currentFig = gcf;
            currentFig.Units = "inches";
            currentFig.Position = [1 1 2.8330 3.6670];
            set(gca,'fontsize',14);
            subplot(2,1,1)
            ylabel('Voltage (mV)')
            plot(1e3*t,1e3*rawSig(:,chan,trial),'linewidth',2,'color','k');
            title('Raw Signal')
            xlim([1000 1015])
            ylim([-0.5 0.5])
            ylabel('Voltage (mV)')
            
            
            subplot(2,1,2)
            currentFig = gcf;
            set(gca,'fontsize',14);
            xlabel('Time (ms)')
            plot(1e3*t',1e3*rawSigTemp,'linewidth',2,'color',[204 85 0]/255);
            title('Recovered Signal')
            xlabel('Time (ms)')
            xlim([1000 1015])
            ylim([-0.5 0.5])
            
            
        end
        
        processedSig(:,chan,trial) = rawSigTemp;
        fprintf(['-------Template Processed - Channel ' num2str(chan) '--' 'Trial ' num2str(trial) '-----\n'])
    end
end
%%
% plot trials belong to particular clusters

if plotIt
    chan = chanInt;
    templateArrayShortened = templateListVec{chan};
    t = 1e3*[0:size(templateArrayShortened,1)-1]/fs;
    figUnsortTime = figure;
    figUnsortTime.Units = "inches";
    figUnsortTime.Position = [1 1 1 10/2];
    numIndices = 20;
    indices = randi(size(templateArrayShortened,2),numIndices,1);
    subselectTrials = templateArrayShortened(:,indices);
    maxSubselect = 1e3*max(abs(subselectTrials(:)));
    for ii = 1:numIndices
        vizFunc.smplot(numIndices,1,ii,'axis','off','right',0.02,'bottom',0.02,'left',0.02,'top',0.02)
        plot(t, 1e3*subselectTrials(:,ii),'color','k','linewidth',2)
        ylim([-maxSubselect maxSubselect])
        set(gca,'XTick',[], 'YTick', [],'YLabel',[], 'Xlabel',[],'Visible','off')
    end
    
    obj = vizFunc.scalebar;
    obj.XLen = 0.5;              %X-Length, 10.
    obj.XUnit = 'ms';            %X-Unit, 'm'.
    obj.YLen = 10;
    obj.YUnit = 'mV';
    set(gca,'fontsize',18)
    %     obj.Position = [1,-1];
    %     obj.hTextX_Pos = [1,-2]; %move only the LABEL position
    %     obj.hTextY_Pos =  [0.5,-1];
    obj.hLineY(2).LineWidth = 5;
    obj.hLineY(1).LineWidth = 5;
    obj.hLineX(2).LineWidth = 5;
    obj.hLineX(1).LineWidth = 5;
    %     obj.Border = 'LL';          %'LL'(default), 'LR', 'UL', 'UR'
    
    templates = templateSubtractCell{chan};
    templatesChoiceTimeSeries = templates(indices);
    figSortTime = figure;
    figSortTime.Units = "inches";
    figSortTime.Position = [1 1 1 10/2];
    tempLabels = labelInfo.tempLabels{chan};
    minLength1 = 7;
    minLength2 = length(unique(tempLabels))+2;
    
    for ii = 1:numIndices
        colorInt = colors(templatesChoiceTimeSeries(ii),:);
        vizFunc.smplot(numIndices,1,ii,'axis','off','right',0.02,'bottom',0.02,'left',0.02,'top',0.02)
        plot(t, 1e3*subselectTrials(:,ii),'color',colorInt,'linewidth',2)
        ylim([-maxSubselect maxSubselect])
        
        set(gca,'XTick',[], 'YTick', [],'YLabel',[], 'Xlabel',[],'Visible','off')
    end
    
    obj = vizFunc.scalebar;
    obj.XLen = 0.5;              %X-Length, 10.
    obj.XUnit = 'ms';            %X-Unit, 'm'.
    obj.YLen = 10;
    obj.YUnit = 'mV';
    set(gca,'fontsize',18)
    %     obj.Position = [1,-1];
    %     obj.hTextX_Pos = [1,-2]; %move only the LABEL position
    %     obj.hTextY_Pos =  [0.5,-1];
    obj.hLineY(2).LineWidth = 5;
    obj.hLineY(1).LineWidth = 5;
    obj.hLineX(2).LineWidth = 5;
    obj.hLineX(1).LineWidth = 5;
    %     obj.Border = 'LL';          %'LL'(default), 'LR', 'UL', 'UR'
    %%
    CTdiv = cbrewerHelper.cbrewer('div', 'BrBG', 50);
    numIndicesISC = 100;
    indicesISC = randi(size(templateArrayShortened,2),numIndicesISC,1);
    templateArrayInt = templateArrayShortened(:,indicesISC);
    
    figHeatMapUnsort = figure;
    figHeatMapUnsort.Units = "inches";
    figHeatMapUnsort.Position = [1 1 2 4];
    imagesc(1e3*templateArrayInt')
    xlabel('Sample')
    ylabel('Trial #')
    set(gca,'fontsize',10)
    colormap(CTdiv)
    caxis([-max(abs(1e3*templateArrayInt(:))) max(abs(1e3*templateArrayInt(:)))])
    c = colorbar();
    c.Label.String = 'Voltage (mV)';
    %     title({'Subselection of Trials',' Before Template Matching'})
    
    %%
    templates = templateSubtractCell{chan};
    templatesChoice = templates(indicesISC);
    figSorted = figure;
    figSorted.Units = "inches";
    figSorted.Position = [1 1 2 4];
    counter = 1;
    for templateInterest = unique(templatesChoice)'
        indicesSelect = find(templatesChoice == templateInterest);
        subplot(length(unique(templatesChoice)),1,counter)
        %         vizFunc.smplot(length(unique(templatesChoice)),1,counter,'axis','off',...
        %             'right',0.1,'bottom',0.1,'left',0.1,'top',0.1)
        
        imagesc(1e3*templateArrayInt(:,indicesSelect)')
        if counter == 1
            %             title({'Sorted Subselection of Trials','By Template'})
        end
        colormap(CTdiv)
        caxis([-max(abs(1e3*templateArrayShortened(:))) max(abs(1e3*templateArrayShortened(:)))])
        if counter < length(unique(templatesChoice))
            set(gca,'YLabel',[], 'Xlabel',[],'Xtick',[])
            
        end
        set(gca,'fontsize',10)
        counter = counter + 1;
        
    end
    xlabel('Sample')
    ylabel('Trial #')
    set(gca,'fontsize',10)
    caxis([-max(abs(1e3*templateArrayShortened(:))) max(abs(1e3*templateArrayShortened(:)))])
    % c = colorbar();
    %c.Label.String = 'Voltage (mV)';
    %%
    templateArrayExtracted = templateArrayCellOutput{chan};
    templateArrayShortened = templateArrayExtracted(maxLocation+bracketRange,:);
    figureTemplate = figure;
    figureTemplate.Units = "inches";
    figureTemplate.Position =[1 1 2 4];
    counter = 1;
    uniqueVec = unique(templatesChoice)';
    t = 1e3*(0:size(templateArrayShortened,1)-1)/fs;
    for templateInterest = uniqueVec(1:end)
        colorInt = colors(counter,:);
        vizFunc.smplot(length(uniqueVec),1,counter,'axis','off','right',0.02,'bottom',0.02,'left',0.02,'top',0.02)
        
        % subplot(length(uniqueVec),1,counter)
        ylim([-maxSubselect maxSubselect])
        
        plot(t,1e3*templateArrayShortened(:,templateInterest),'color',colorInt,'linewidth',2)
        counter = counter + 1;
        set(gca,'fontsize',14)
        set(gca,'XTick',[], 'YTick', [],'YLabel',[], 'Xlabel',[],'Visible','off')
        
    end
    xlabel('Time (ms')
    ylabel('Voltage (ms)')
    
    obj = vizFunc.scalebar;
    obj.XLen = 0.5;              %X-Length, 10.
    obj.XUnit = 'ms';            %X-Unit, 'm'.
    obj.YLen = 10;
    obj.YUnit = 'mV';
    set(gca,'fontsize',18)
    %     obj.Position = [1,-1];
    %     obj.hTextX_Pos = [1,-2]; %move only the LABEL position
    %     obj.hTextY_Pos =  [0.5,-1];
    obj.hLineY(2).LineWidth = 5;
    obj.hLineY(1).LineWidth = 5;
    obj.hLineX(2).LineWidth = 5;
    obj.hLineX(1).LineWidth = 5;
    
    figure
    clusterer = clusterInt;
    currentFig = gcf;
    currentFig.Units = "inches";
    currentFig.Position = [1 1 1.8 1.8];
    h = scatter(1e3*clusterer.data(:,maxSub),1e3*clusterer.data(:,maxSub+1),20,'filled');
    %ylabel('time point 1 : voltage (mV)')
    ylabel('V(t=2)')
    xlabel('V(t=1)')
    % xlabel('time point 2 : voltage (mV)')
    set( h.Parent,'tickdir','out','box','off' );
    
    tempLabels = clusterer.labels;
    tempLabels(tempLabels == 0) = max(unique(tempLabels))+1;
    [C,ia,ic] = unique(tempLabels);
    labelInfo.C{chan} = C;
    labelInfo.ic{chan} = ia;
    labelInfo.tempLabels{chan} = tempLabels;
    h.CData = ic;
    colormap(colorsNew);
    labelInfo.colors{chan} = colors;
    set(gca,'fontsize',12);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end