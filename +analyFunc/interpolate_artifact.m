function [processedSig] = interpolate_artifact(rawSig,varargin)
%USAGE:
% This function will perform an interpolation scheme for artifacts on a
% trial by trial, channel by channel basis, implementing either a linear
% interpolation scheme, or a pchip interpolation scheme
%
% raw_sig = samples x channels x trials
% pre = the number of ms before which the stimulation pulse onset as
% detected by a thresholding method should still be considered as artifact
% post = the number of ms after which the stimulation pulse onset as
% detected by a thresholding method should still be considered as artifact
% pre_interp = the number of ms before the stimulation which to consider an
% interpolation scheme on. Does not apply to the linear case
% post_interp = the number of ms before the stimulation which to consider an
% interpolation scheme on. Does not apply to the linear case
% fs = sampling rate (Hz)
% plotIt = plot intermediate steps if true

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;

validData = @(x) isnumeric(x) && size(x,3)>2;
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
postInterp = p.Results.postInterp;
stimChans = p.Results.stimChans;
bads = p.Results.bads;
fixedDistance = p.Results.fixedDistance;
fs = p.Results.fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define matrix of zeros
processedSig = zeros(size(rawSig));

% make a vector of the good channels to process
numChans = size(rawSig,2);
[goods,goodVec] = helpFunc.goodChannel_extract('bads',bads,'stimchans',stimChans,...,
    'numChans',numChans);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(['-------Interpolation-------- \n'])
fprintf(['-------' type '-------- \n'])

[startInds,endInds] = analyFunc.get_artifact_indices(rawSig,'pre',pre,'post',post,'plotIt',...,
    plotIt,'fixedDistance',fixedDistance,'fs',fs,'goodVec',goodVec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for trial = 1:size(rawSig,3)
    
    for chan = goodVec
        rawSigTemp = rawSig(:,chan,trial);
        for sts = 1:length(startInds{trial})
            win = startInds{trial}(sts):endInds{trial}(sts);
            switch type
                case 'linear'
                    rawSigTemp(win) = interp1([startInds{trial}(sts)-1 endInds{trial}(sts)+1],...
                        rawSigTemp([startInds{trial}(sts)-1 endInds{trial}(sts)+1]), startInds{trial}(sts):endInds{trial}(sts));
                case 'pchip'
                    preInterpSamps = round(preInterp*fs/1e3);
                    postInterpSamps = round(postInterp*fs/1e3);
                    rawSigTemp(win) = interp1([startInds{trial}(sts)-preInterpSamps:startInds{trial}(sts)-1 endInds{trial}(sts):endInds{trial}(sts)+postInterpSamps],...
                        rawSigTemp([startInds{trial}(sts)-preInterpSamps:startInds{trial}(sts)-1 endInds{trial}(sts):endInds{trial}(sts)+postInterpSamps]),...
                        startInds{trial}(sts):endInds{trial}(sts),'pchip');
            end
            
        end
        
        if plotIt && (trial == 10 || trial == 1000)
            figure
            plot(rawSigTemp,'linewidth',2)
            hold on
            plot(rawSig(:,chan,trial),'linewidth',2)
            vline(startInds{trial})
            vline(endInds{trial},'g')
        end
        
        processedSig(:,chan,trial) = rawSigTemp;
    end
end

fprintf(['-------Finished-------- \n \n'])

end