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
addParameter(p,'amntPreAverage',3,@isnumeric);
addParameter(p,'bracketRange',[-8:8],@isnumeric);

addParameter(p,'expThreshVoltageCut',95,@isnumeric);
addParameter(p,'expThreshDiffCut',95,@isnumeric);

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
amntPreAverage = p.Results.amntPreAverage;
maxLocation = p.Results.maxLocation;
bracketRange = p.Results.bracketRange;

expThreshVoltageCut = p.Results.expThreshVoltageCut;
expThreshDiffCut = p.Results.expThreshDiffCut;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
templateArrayCellOutput = {};
processedSig = zeros(size(rawSig));
templateSubtractCell = {};
templateListVec = {};

fprintf(['-------Dictionary-------- \n'])

plotIt = 1;

for chan = goodVec
    templateArray = templateArrayCell{chan};
    
    % extract max amplitude for a given channel
    maxAmpsChan = max(maxAmps(chan,:));
    
    % shorten data to be centered around the peak +/- the bracketRange. In
    % this way there is less clustering around non-discriminative data
    % points.
    templateArrayShortened = templateArray(maxLocation+bracketRange,:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % data is in "templateArrayShorted". We will initiate a new HDBSCAN instance
    clusterer = HDBSCAN.HDBSCAN( templateArrayShortened');
    
    % we can view our data matrix size
    fprintf( 'Number of points: %i\n',clusterer.nPoints );
    fprintf( 'Number of dimensions: %i\n',clusterer.nDims );
    
    try
        % (1) directly set the parameters
        clusterer.minpts = 2;
        clusterer.minclustsize = 3;
        clusterer.outlierThresh = 0.90;
        
        clusterer.metric = distanceMetricDbscan;
        clusterer.fit_model(); 			% trains a cluster hierarchy
    catch
        % (1) directly set the parameters
        clusterer.minpts = 3;
        clusterer.minclustsize = 3;
        clusterer.outlierThresh = 0.90;
        
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
    
    if plotIt && chan == 28
        CT = cbrewerHelper.cbrewer('qual', 'Dark2', 8);
        
        repeats = max( 1,ceil( max( clusterer.labels )/8 ) );
        colors = [CT(1,:);repmat( CT(2:end,:),repeats,1 )];
        
        figure
        h = scatter(1e3*clusterer.data(:,maxLocation),1e3*clusterer.data(:,maxLocation-1),'.');
        ylabel('time point 1 : voltage (mV)')
        xlabel('time point 2 : voltage (mV)')
        set( h.Parent,'tickdir','out','box','off' );
        h.CData = clusterer.labels;
        colormap(colors);
        
        set(gca,'fontsize',18)
        
        figure
        subplot(2,1,1)
        h = scatter(1e3*clusterer.data(:,maxLocation),1e3*clusterer.data(:,maxLocation-1),'.');
        ylabel('time point 1 : voltage (mV)')
        xlabel('time point 2 : voltage (mV)')
        set( h.Parent,'tickdir','out','box','off' );
        h.CData = clusterer.labels;
        colormap(colors);
        
        set(gca,'fontsize',18)
        
        subplot(2,1,2)
        labelsInt = unique(clusterer.labels);
        labelsInt = labelsInt(labelsInt~=0);
        t = 1e3*[0:size(templateArrayShortened,1)-1]/fs;
        plot(t,1e3*templateArrayShortened,'color',[0.5 0.5 0.5])
        hold on
        h2 = plot(t,1e3*templateArrayExtracted(maxLocation+bracketRange,:),'linewidth',2);
        set(h2, {'color'}, num2cell(colors(labelsInt,:),2));
        
        colormap( colors );
        
        ylabel('Voltage (mV)')
        xlabel('Time (ms)')
        set(gca,'fontsize',18)
        title('Dictionary of Templates from Raw Traces')
        
        index = 1;
        for labelsIndiv = labelsInt'
            figure
            t = 1e3*[0:size(templateArrayShortened,1)-1]/fs;
            plot(t,1e3*templateArrayShortened(:,clusterer.labels==labelsIndiv),'color',[0.5 0.5 0.5])
            hold on
            h2 = plot(t,1e3*templateArrayExtracted(maxLocation+bracketRange,index),'linewidth',2);
            set(h2, {'color'}, num2cell(colors(labelsIndiv,:),2));
            ylabel('Voltage (mV)')
            xlabel('Time (ms)')
            set(gca,'fontsize',18)
            title('Dictionary of Templates from Raw Traces')
            index = index + 1;
        end
        
        templateListVec{chan} = templateArrayShortened;
        
    end
    
    % assign templates to channel
    templateArrayCellOutput{chan} = templateArrayExtracted;
    fprintf(['-------Artifact Channel ' num2str(chan) ' -------- \n'])
    
end

fprintf(['-------Finished clustering artifacts-------- \n'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% now do the template subtraction
for trial = 1:size(rawSig,3)
    
    for chan = goodVec
        
        rawSigTemp = rawSig(:,chan,trial);
        templates = templateArrayCellOutput{chan};
        
        % add on the trial one
        templates = [templates mean(templateTrial{chan}{trial},2)];
        
        % ensure no subtraction of exponential
        if recoverExp
            templates = analyFunc.recover_EP(templates,fs,'threshDiffCut',expThreshDiffCut,'threshVoltageCut',expThreshVoltageCut);
        end
        
        for sts = 1:length(startInds{trial}{chan})
            win = startInds{trial}{chan}(sts):endInds{trial}{chan}(sts);
            extractedSig = rawSigTemp(win);
            extractedSig = extractedSig - mean(extractedSig(1:amntPreAverage));
            
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
            
            templatesStsShortened = templatesSts(maxLocation+bracketRange,:);
            extractedSigShortened = extractedSig(maxLocation+bracketRange,:);
            
            switch distanceMetricSigMatch
                case 'correlation'
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
            % which template won out
            if isempty(templateSubtractCell{chan}
                templateSubtractCell{chan} = index;
            else
                templateSubtractCell{chan} = [templateSubtractCell{chan}; index];
            end
            rawSigTemp(win) = rawSigTemp(win) - templateSubtract;
            
            if plotIt
                if 1 && chan == 28 && (sts == 1 || sts == 2 || sts == 10)
                    figure
                    t = [0:length(extractedSig)-1]/fs;
                    
                    plot(1e3*t,1e3*extractedSig,'linewidth',2,'color','k')
                    hold on
                    plot(1e3*t,1e3*templateSubtract,'linewidth',2,'color',[0.5 0.5 0.5])
                    plot(1e3*t,1e3*extractedSig-templateSubtract,'linewidth',2,'color',[204 85 0]/255)
                    legend('Raw Signal','Template Selected','Recovered Signal');
                    set(gca,'fontsize',18)
                    title('Raw signal, Template, and Recovered signal')
                    ylabel('Voltage (mV)')
                    xlabel('Time (ms)')
                end
            end
        end
        
        if plotIt && (trial == 10 || trial == 15 || trial == 20) && chan == 28
            figure
            t = [0:length(rawSigTemp)-1]/fs;
            plot(1e3*t',1e3*rawSigTemp,'linewidth',2,'color',[204 85 0]/255)
            hold on
            plot(1e3*t,1e3*rawSig(:,chan,trial),'linewidth',2,'color','k')
            for indexPlot = 1:length(startInds{trial}{chan})
                tempBox = vizFunc.highlight(gca, [1e3*startInds{trial}{chan}(indexPlot)/fs 1e3*endInds{trial}{chan}(indexPlot)/fs], [1e3*min(rawSig(:,chan,trial))-10 1e3*max(rawSig(:,chan,trial))+10], [.5 .5 .5]);
            end
            
            legend({'Recovered Signal','Raw Signal','Artifact Windows'})
            set(gca,'fontsize',18)
            title('Raw vs. Recovered Sig')
            xlabel('Time (ms)')
            ylabel('Voltage (mV)')
        end
        
        processedSig(:,chan,trial) = rawSigTemp;
        fprintf(['-------Template Processed - Channel ' num2str(chan) '--' 'Trial ' num2str(trial) '-----\n'])
    end
    
end

% plot trials belong to particular clusters

if plotIt && chan == 28
    
    templateArrayShortened = templateListVec{chan};
    t = 1e3*[0:size(templateArrayShortened,1)-1]/fs;
    figure
    numIndices = 50;
    indices = randi(size(templateArrayShortened,2),numIndices,1);
    for ii = 1:numIndices
        vizFunc.smplot(32,1,ii,'axis','off')
        plot(t,templateArrayShortened(:,indices(ii)),'color','k','linewidth',2)
        set(gca,'XTick',[], 'YTick', [],'YLabel',[], 'Xlabel',[],'Visible','off')
    end
    
    obj = vizFunc.scalebar;
    obj.XLen = 200;              %X-Length, 10.
    obj.XUnit = 'ms';            %X-Unit, 'm'.
    obj.YLen = 3;
    obj.YUnit = 'mV';
    
    obj.Position = [20,-130];
    obj.hTextX_Pos = [5,-50]; %move only the LABEL position
    obj.hTextY_Pos =  [-45,-40];
    obj.hLineY(2).LineWidth = 5;
    obj.hLineY(1).LineWidth = 5;
    obj.hLineX(2).LineWidth = 5;
    obj.hLineX(1).LineWidth = 5;
    obj.Border = 'LL';          %'LL'(default), 'LR', 'UL', 'UR'
    
    CTdiv = cbrewerHelper.cbrewer('div', 'BrBG', 50);
    numIndicesISC = 100;
    indicesISC = randi(size(templateArrayShortened,2),numIndicesISC,1);
    figure
    imagesc(1e3*templateArrayShortened(:,indicesISC)')
    xlabel('Time (ms)')
    ylabel('Trial #')
    set(gca,'fontsize',18)
    colormap(CTdiv)
    caxis([-max(abs(1e3*templateArrayShortened(:))) max(abs(1e3*templateArrayShortened(:)))])
    c = colorbar();
    c.Label.String = 'Voltage (mV)';
    
    
    templates = templateArrayCellOutput{chan};
    figure
    imagesc(1e3*templateArrayShortened(:,indicesISC)')
    xlabel('Time (ms)')
    ylabel('Trial #')
    set(gca,'fontsize',18)
    colormap(CTdiv)
    caxis([-max(abs(1e3*templateArrayShortened(:))) max(abs(1e3*templateArrayShortened(:)))])
    c = colorbar();
    c.Label.String = 'Voltage (mV)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end