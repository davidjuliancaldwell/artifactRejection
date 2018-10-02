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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
templateArrayCellOutput = {};
processedSig = zeros(size(rawSig));


fprintf(['-------Dictionary-------- \n'])


for chan = goodVec
    templateArray = templateArrayCell{chan};
    
    % extract max amplitude for a given channel
    maxAmpsChan = max(maxAmps(chan,:));
    
    templateArrayShortened = templateArray(maxLocation+bracketRange,:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Assume our data is in the n x m matrix "X". We will initiate a new HDBSCAN instance
    clusterer = HDBSCAN.HDBSCAN( templateArrayShortened');
    
    % we can view our data matrix size
    fprintf( 'Number of points: %i\n',clusterer.nPoints );
    fprintf( 'Number of dimensions: %i\n',clusterer.nDims );
    
    % (1) directly set the parameters
    clusterer.minpts = 3;
    clusterer.minclustsize = 5;
    clusterer.outlierThresh = 0.90;
    clusterer.metric = 'eucl';
    clusterer.fit_model(); 			% trains a cluster hierarchy
    clusterer.get_best_clusters(); 	% finds the optimal "flat" clustering scheme
    clusterer.get_membership();		% assigns cluster labels to the points in X
    
    % (2) call run_hdbscan() with optional inputs. This is the prefered/easier method
    %clusterer.run_hdbscan( 10,20,[],0.85 );
    
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
    
    
    if plotIt
        figure
        clusterer.plot_tree();
        
        % we can also visualize the actual clusters in a 2D or 3D space (depending on self.nDims)
        figure
        subplot(3,1,1)
        clusterer.plot_clusters([6,7,8]); % plots a scatter of the points, color coded by the associated labels
        %figure
        subplot(3,1,2)
        clusterer.plot_clusters( [9,10,11] ); % specifies the scatter to use the 1st, 4th, and 5th columns of X
        
        subplot(3,1,3)
        
        plot(templateArrayExtracted)
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
        % templates = [templates mean(templateTrial{chan}{trial},2)];
        
        % ensure no subtraction of exponential
        if recoverExp
            templates = analyFunc.recover_EP(templates,fs);
        end
        
        for sts = 1:length(startInds{trial}{chan})
            win = startInds{trial}{chan}(sts):endInds{trial}{chan}(sts);
            extractedSig = rawSigTemp(win);
            extractedSig = extractedSig - mean(extractedSig(1:amntPreAverage));
            
            % find best artifact
            % get them to be the same length
            templatesSts = templates(1:length(extractedSig),:);
            
            
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
            
            rawSigTemp(win) = rawSigTemp(win) - templateSubtract;
            if plotIt
                if 1 && chan == 28 && (sts == 1 || sts == 2 || sts == 3)
                    figure
                    plot(extractedSig)
                    hold on
                    plot(templateSubtract)
                    plot(extractedSig-templateSubtract)
                    legend('extracted','template','subtracted');
                end
            end
        end
        
        if plotIt
            if 1 && chan == 28
                figure
                plot(rawSigTemp,'linewidth',2)
                hold on
                plot(rawSig(:,chan,trial),'linewidth',2)
                
                vline(startInds{trial}{chan})
                vline(endInds{trial}{chan},'g')
                xlim([1.221e4 1.236e4])
                
            end
        end
        processedSig(:,chan,trial) = rawSigTemp;
        fprintf(['-------Template Processed - Channel ' num2str(chan) '--' 'Trial ' num2str(trial) '-----\n'])
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end