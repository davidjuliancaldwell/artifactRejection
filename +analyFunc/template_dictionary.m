function [processedSig,templateArrayCellOutput] = template_dictionary(templateArrayCell,templateTrial,rawSig,varargin)
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
%  preInterp - the number of ms before the stimulation which to consider an
%              interpolation scheme on. Does not apply to the linear case
% postInterp - the number of ms before the stimulation which to consider an
%              interpolation scheme on. Does not apply to the linear case
%          fs - sampling rate (Hz)
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

addParameter(p,'plotIt',0,@(x) x==0 || x ==1);

addParameter(p,'distanceMetricDbscan','eucl',@isstr);
addParameter(p,'distanceMetricSigMatch','eucl',@isstr);

addParameter(p,'goodVec',[1:64],@isnumeric);
addParameter(p,'startInds',[],@iscell);
addParameter(p,'endInds',[],@iscell);

addParameter(p,'recoverExp',1,@(x) x==0 || x ==1);


p.parse(templateArrayCell,templateTrial,rawSig,varargin{:});

templateArrayCell = p.Results.templateArrayCell;
templateTrial = p.Results.templateTrial;
rawSig = p.Results.rawSig;

plotIt = p.Results.plotIt;
distanceMetricDbscan = p.Results.distanceMetricDbscan;
distanceMetricSigMatch = p.Results.distanceMetricSigMatch;
goodVec = p.Results.goodVec;
startInds = p.Results.startInds;
endInds = p.Results.endInds;

recoverExp = p.Results.recoverExp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

templateArrayCellOutput = {};
processedSig = zeros(size(rawSig));

fprintf(['-------Dictionary-------- \n'])

for chan = goodVec
    templateArrayExtracted = [];
    templateArray = templateArrayCell{chan};
    if strcmp(distanceMetricDbscan,'eucl')
        if max(templateArray(:)) < 3e-4
            distanceDBscan = 1e-3;
        else
            distanceDBscan = 1e-4;
        end
    else
        if max(templateArray(:)) < 3e-4
            distanceDBscan = 0.9;
        else
            distanceDBscan = 0.95;
        end
    end
    
    [c,ptsC,centres] = analyFunc.db_scan(templateArray,distanceDBscan,1,distanceMetricDbscan);
    
    vectorUniq = unique(ptsC);
    if plotIt
        figure
        hold on
        vectorUniq = unique(ptsC)';
        for i = 1:length(vectorUniq)
            plot(mean(templateArray(:,ptsC==i),2));
        end
    end
    
    
    for i = 1:length(vectorUniq)
        meanTempArray = mean(templateArray(:,ptsC==i),2);
        templateArrayExtracted = [templateArrayExtracted (meanTempArray )]; %no subtraction
    end
    
    templateArrayCellOutput{chan} = templateArrayExtracted;
    
    fprintf(['-------Artifact Channel ' num2str(chan) ' -------- \n'])
    
    
end

fprintf(['-------Finished clustering artifacts-------- \n'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% now do the template subtraction
for trial = 1:size(rawSig,3)
    
    for chan = goodVec
        
        rawSigTemp = rawSig(:,chan,trial);
        templates = templateArrayCellOutput{chan};
        
        % add on the trial one
        templates = [templates mean(templateTrial{chan}{trial},2)];
        
        
        % ensure no subtraction of exponential
        
        if recoverExp
            templates = analyFunc.recover_EP(templates);
        end
        
        for sts = 1:length(startInds{trial}{chan})
            win = startInds{trial}{chan}(sts):endInds{trial}{chan}(sts);
            extractedSig = rawSigTemp(win);
            extractedSig = extractedSig - extractedSig(1);
            
            % find best artifact
            % get them to be the same length
            templatesSts = templates(1:length(extractedSig),:);
            
            switch distanceMetricSigMatch
                case 'correlation'
                    % correlation
                    correlations = corr(extractedSig,templatesSts);
                    [~,index] = max((correlations));
                    
                case 'cosineSim'
                    % - cosine similarity
                    denominator = sqrt(sum(extractedSig.*extractedSig)).*sqrt(sum(templatesSts.*templates));
                    numerator = (extractedSig'*templatesSts);
                    correlations = numerator./denominator;
                    [~,index] = max(abs(correlations));
                    
                case 'eucl'
                    % distance
                    v = templatesSts - repmat(extractedSig,1,size(templatesSts,2));
                    distance = sum(v.*v);
                    [~,index] = min(distance);
            end
            
            templateSubtract = templatesSts(:,index);
            
            %             k = 2;
            %             all = templateArrayCell{chan};
            %             [u,s,v] = svd(all);
            %
            %             c = u(:,1:k);
            %             %c = template_subtract;
            %             a = rawSigTemp(win);
            %             d = (c'*c)\(c'*a);
            %             clean = a - (c'*c)\(c'*a);
            %
            rawSigTemp(win) = rawSigTemp(win) - templateSubtract;
        end
        
        if plotIt && (trial == 10 || trial == 1000)
            figure
            plot(rawSigTemp,'linewidth',2)
            hold on
            plot(rawSig(:,chan,trial),'linewidth',2)
            
            vline(startInds{trial}{chan})
            vline(endInds{trial}{chan},'g')
            
            figure
            plot(extractedSig)
            hold on
            plot(templateSubtract)
            plot(extractedSig-templateSubtract)
            legend('extracted','template','subtracted');
        end
        processedSig(:,chan,trial) = rawSigTemp;
            fprintf(['-------Template Processed - Channel ' num2str(chan) '--' 'Trial ' num2str(trial) '-----\n'])
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end