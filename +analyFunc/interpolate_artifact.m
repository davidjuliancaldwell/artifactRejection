function [processedSig,startInds,endInds] = interpolate_artifact(rawSig,varargin)
% Usage:  [processedSig,startInds,endInds] = interpolate_artifact(rawSig,varargin)
%
% This function will perform an interpolation scheme for artifacts on a
% trial by trial, channel by channel basis, implementing either a linear
% interpolation scheme, or a pchip interpolation scheme
%
% Arguments:
%   Required:
%   rawSig - samples x channels x trials
%
%   Optional:
%        pre - The number of ms before which the stimulation pulse onset as
%              detected by a thresholding method should still be considered
%              as artifact
%       post - The number of ms after which the stimulation pulse onset as
%              detected by a thresholding method should still be
%              considered as artifact
%  preInterp - The number of ms before the stimulation which to consider an
%              interpolation scheme on. Does not apply to the linear case
% postInterp - The number of ms before the stimulation which to consider an
%              interpolation scheme on. Does not apply to the linear case
%          fs - Sampling rate (Hz)
%      plotIt - Plot intermediate steps if true
% useFixedEnd - Use a fixed end
%
% Returns:
%   processedSig - Processed signal following interpolation scheme
%      startInds - Cell array of the start indices each artifact for each
%      channel and trial - startInds{trial}{channel}
%       endsInds - Cell array of the end indices of each artifact for each
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

validData = @(x) isnumeric(x);
addRequired(p,'rawSig',validData);

addParameter(p,'plotIt',0,@(x) x==0 || x ==1);
addParameter(p,'useFixedEnd',0,@(x) x==0 || x ==1);

addParameter(p,'type','linear',@isstr);
addParameter(p,'pre',0.4096,@isnumeric);
addParameter(p,'post',0.4096,@isnumeric);
addParameter(p,'preInterp',0.2,@isnumeric);
addParameter(p,'postInterp',0.2,@isnumeric);
addParameter(p,'stimChans',[],@isnumeric);
addParameter(p,'bads',[],@isnumeric);
addParameter(p,'fixedDistance',2,@isnumeric);
addParameter(p,'fs',12207,@isnumeric);

p.parse(rawSig,varargin{:});
rawSig = p.Results.rawSig;
plotIt = p.Results.plotIt;
type = p.Results.type;
useFixedEnd = p.Results.useFixedEnd;
pre = p.Results.pre;
post = p.Results.post;
preInterp = p.Results.preInterp;
postInterp = p.Results.postInterp;close all

stimChans = p.Results.stimChans;
bads = p.Results.bads;
fixedDistance = p.Results.fixedDistance;
fs = p.Results.fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define matrix of zeros
processedSig = zeros(size(rawSig));

% make a vector of the good channels to process
numChans = size(rawSig,2);
[goods,goodVec] = helpFunc.good_channel_extract('bads',bads,'stimchans',stimChans,...
    'numChans',numChans);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(['-------Interpolation-------- \n'])
fprintf(['-------' type '-------- \n'])

[startInds,endInds] = analyFunc.get_artifact_indices(rawSig,'pre',pre,'post',post,'plotIt',...,
    0,'useFixedEnd',useFixedEnd,'fixedDistance',fixedDistance,'fs',fs,'goodVec',goodVec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for trial = 1:size(rawSig,3)
    
    for chan = goodVec
        rawSigTemp = rawSig(:,chan,trial);
        
        for sts = 1:length(startInds{trial}{chan})
            
            win = startInds{trial}{chan}(sts):endInds{trial}{chan}(sts);
            switch type
                case 'linear'
                    %%
                    rawSigTemp(win) = interp1([startInds{trial}{chan}(sts)-1 endInds{trial}{chan}(sts)+1],...
                        rawSigTemp([startInds{trial}{chan}(sts)-1 endInds{trial}{chan}(sts)+1]), startInds{trial}{chan}(sts):endInds{trial}{chan}(sts));
                case 'pchip'
                    preInterpSamps = round(preInterp*fs/1e3);
                    postInterpSamps = round(postInterp*fs/1e3);
                    rawSigTemp(win) = interp1([startInds{trial}{chan}(sts)-preInterpSamps:startInds{trial}{chan}(sts)-1 endInds{trial}{chan}(sts):endInds{trial}{chan}(sts)+postInterpSamps],...
                        rawSigTemp([startInds{trial}{chan}(sts)-preInterpSamps:startInds{trial}{chan}(sts)-1 endInds{trial}{chan}(sts):endInds{trial}{chan}(sts)+postInterpSamps]),...
                        startInds{trial}{chan}(sts):endInds{trial}{chan}(sts),'pchip');
            end
        end
        
        if plotIt && (trial == 10 || trial == 1000)
            figure
            t = [0:length(rawSigTemp)-1]/fs;
            plot(1e3*t',1e6*rawSigTemp,'linewidth',2,'color','r')
            hold on
            plot(1e3*t,1e6*rawSig(:,chan,trial),'linewidth',2,'color','k')
            for indexPlot = 1:length(startInds{trial}{chan})
               tempBox = vizFunc.highlight(gca, [1e3*startInds{trial}{chan}(indexPlot)/fs 1e3*endInds{trial}{chan}(indexPlot)/fs], [1e6*min(rawSig(:,chan,trial)) 1e6*max(rawSig(:,chan,trial))], [.5 .5 .5]);
            end
     
            legend({'linear interpolation','raw signal','artifact window'})
            set(gca,'fontsize',18)
            title('Raw vs. Processed Sig')
            xlabel('Time (ms)')
            ylabel('Voltage (\muV)')
        end
        
        processedSig(:,chan,trial) = rawSigTemp;
        fprintf(['-------Interpolation - Channel ' num2str(chan) '--' 'Trial ' num2str(trial) '-----\n'])
        
    end
end

fprintf(['-------Finished-------- \n \n'])

end