function [processedSig,templateArrayCellOutput,templateTrial,startInds,endInds] = template_subtract(rawSig,varargin)
%USAGE:
% This function will perform a template subtraction scheme for artifacts on
% a trial by trial, channel by channel basis. This function will build up
% a dictionary of artifacts, best match the template, and
%
% raw_sig = samples x channels x trials
% pre = the number of ms before which the stimulation pulse onset as
% detected by a thresholding method should still be considered as artifact
% post = the number of ms after which the stimulation pulse onset as
% detected by a thresholding method should still be considered as artifact
% fs = sampling rate (Hz)
% plotIt = plot intermediate steps if true

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get inputs
p = inputParser;

validData = @(x) isnumeric(x) && size(x,3)>2;
addRequired(p,'rawSig',validData);

addParameter(p,'plotIt',0,@(x) x==0 || x ==1);
addParameter(p,'useFixedEnd',0,@(x) x==0 || x ==1);

addParameter(p,'type','average',@isstr);
addParameter(p,'distanceMetricDbscan','eucl',@isstr);
addParameter(p,'distanceMetricSigMatch','eucl',@isstr);

addParameter(p,'pre',0.4096,@isnumeric);
addParameter(p,'post',0.4096,@isnumeric);
addParameter(p,'preInterp',0.2,@isnumeric);
addParameter(p,'postInterp',0.2,@isnumeric);
addParameter(p,'stimChans',[],@isnumeric);
addParameter(p,'bads',[],@isnumeric);
addParameter(p,'fixedDistance',2,@isnumeric);
addParameter(p,'fs',12207,@isnumeric);
addParameter(p,'amntPreAverage',5,@isnumeric);
addParameter(p,'normalize','firstSamp',@isstr);
addParameter(p,'recoverExp',1,@(x) x==0 || x ==1);

p.parse(rawSig,varargin{:});

rawSig = p.Results.rawSig;
plotIt = p.Results.plotIt;
useFixedEnd = p.Results.useFixedEnd;

type = p.Results.type;
distanceMetricDbscan = p.Results.distanceMetricDbscan;
distanceMetricSigMatch = p.Results.distanceMetricSigMatch;

pre = p.Results.pre;
post = p.Results.post;
preInterp = p.Results.preInterp;
postInterp = p.Results.postInterp;
stimChans = p.Results.stimChans;
bads = p.Results.bads;
fixedDistance = p.Results.fixedDistance;
fs = p.Results.fs;

amntPreAverage = p.Results.amntPreAverage;
normalize = p.Results.normalize;

recoverExp = p.Results.recoverExp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define matrix of zeros
processedSig = zeros(size(rawSig));

% make a vector of the good channels to process
numChans = size(rawSig,2);
[goods,goodVec] = helpFunc.good_channel_extract('numChans',numChans,'bads',bads,'stimChans',stimChans);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get beginnings and ends of artifacts

[startInds,endInds] = analyFunc.get_artifact_indices(rawSig,'pre',pre,'post',post,'plotIt',...,
    plotIt,'useFixedEnd',useFixedEnd,'fixedDistance',fixedDistance,'fs',fs,'goodVec',goodVec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% extract artifacts

[lengthMax,templateCell] = analyFunc.get_artifacts(rawSig,'goodVec',goodVec,...,
    'startInds',startInds,'endInds',endInds,'plotIt',plotIt,'normalize',normalize,'amntPreAverage',amntPreAverage);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get templates all same length

[templateTrial,templateArrayCell] = analyFunc.template_equalize_length(templateCell,rawSig,'lengthMax',...,
    lengthMax,'startInds',startInds,'goodVec',goodVec);

%% build up dictionary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch type
    case 'dictionary'
        [processedSig,templateArrayCellOutput] = analyFunc.template_dictionary(templateArrayCell,templateTrial,rawSig,'plotIt',plotIt,'distanceMetricDbscan',distanceMetricDbscan,...,
            'distanceMetricSigMatch',distanceMetricSigMatch,'goodVec',goodVec,'startInds',startInds,'endInds',endInds,'recoverExp',recoverExp);
        
    case 'average'
        [processedSig,templateArrayCellOutput] = analyFunc.template_average(templateArrayCell,rawSig,'plotIt',plotIt...,
            ,'goodVec',goodVec,'startInds',startInds,'endInds',endInds);
    case 'trial'
        [processedSig,templateArrayCellOutput] = analyFunc.template_trial(templateTrial,rawSig,'plotIt',plotIt...,
            ,'goodVec',goodVec,'startInds',startInds,'endInds',endInds);
end

fprintf(['-------Extracting data-------- \n \n'])

end