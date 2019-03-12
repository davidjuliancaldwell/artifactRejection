function [processedSig,templateArrayCellOutput,templateTrial,startInds,endInds] = template_subtract(rawSig,varargin)
%USAGE:
% This function will perform a template subtraction scheme for artifacts on
% a trial by trial, channel by channel basis. This function will build up
% a dictionary of artifacts, best match the template, and
%
% rawSig = samples x channels x trials
% pre = the number of ms before which the stimulation pulse onset as
% detected by a thresholding method should still be considered as artifact
% post = the number of ms after which the stimulation pulse onset as
% detected by a thresholding method should still be considered as artifact
% fs = sampling rate (Hz)
% plotIt = plot intermediate steps if true

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get inputs
p = inputParser;

validData = @(x) isnumeric(x);
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
addParameter(p,'recoverExp',0,@(x) x==0 || x ==1);
addParameter(p,'minDuration',0,@isnumeric);
addParameter(p,'bracketRange',[-8:8],@isnumeric);
addParameter(p,'onsetThreshold',1.5,@isnumeric);

addParameter(p,'threshVoltageCut',75,@isnumeric);
addParameter(p,'threshDiffCut',75,@isnumeric);

addParameter(p,'expThreshVoltageCut',75,@isnumeric);
addParameter(p,'expThreshDiffCut',75,@isnumeric);

addParameter(p,'chanInt',1,@isnumeric);


addParameter(p,'minPts',2,@isnumeric);
addParameter(p,'minClustSize',3,@isnumeric);
addParameter(p,'outlierThresh',0.95,@isnumeric);

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

onsetThreshold = p.Results.onsetThreshold;
amntPreAverage = p.Results.amntPreAverage;
normalize = p.Results.normalize;

recoverExp = p.Results.recoverExp;
minDuration = p.Results.minDuration;

bracketRange = p.Results.bracketRange; 

threshVoltageCut = p.Results.threshVoltageCut;
threshDiffCut = p.Results.threshDiffCut;

expThreshVoltageCut = p.Results.expThreshVoltageCut;
expThreshDiffCut = p.Results.expThreshDiffCut;

chanInt = p.Results.chanInt;

minPts = p.Results.minPts;
minClustSize = p.Results.minClustSize;
outlierThresh = p.Results.outlierThresh;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define matrix of zeros
processedSig = zeros(size(rawSig));

% make a vector of the good channels to process
numChans = size(rawSig,2);
[goods,goodVec] = helpFunc.good_channel_extract('numChans',numChans,'bads',bads,'stimChans',stimChans);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get beginnings and ends of artifacts

[startInds,endInds] = analyFunc.get_artifact_indices(rawSig,'pre',pre,'post',post,'plotIt',...,
    plotIt,'useFixedEnd',useFixedEnd,'fixedDistance',fixedDistance,'fs',fs,'goodVec',goodVec,...
    'minDuration',minDuration,'threshVoltageCut',threshVoltageCut,'threshDiffCut',threshDiffCut,'onsetThreshold',onsetThreshold);

      %       corrected_signal = artifact_correction(squeeze(processedSig(:,:,1)), index, stimulation_mode, stimulation_length, sampling_frequency, default, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% extract artifacts

[templateCell,lengthMax,maxAmps,maxLocation] = analyFunc.get_artifacts(rawSig,'goodVec',goodVec,...,
    'startInds',startInds,'endInds',endInds,'plotIt',plotIt,'normalize',normalize,'amntPreAverage',amntPreAverage);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get templates all same length

[templateTrial,templateArrayCell] = analyFunc.template_equalize_length(templateCell,rawSig,'lengthMax',...,
    lengthMax,'startInds',startInds,'goodVec',goodVec);

%% build up dictionary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch type
    case 'dictionary'
        [processedSig,templateArrayCellOutput] = analyFunc.template_dictionary(templateArrayCell,templateTrial,rawSig,fs,'plotIt',plotIt,...
            'distanceMetricDbscan',distanceMetricDbscan,'distanceMetricSigMatch',distanceMetricSigMatch,...
            'goodVec',goodVec,'startInds',startInds,'endInds',endInds,'recoverExp',recoverExp,'maxAmps',maxAmps,...
            'amntPreAverage',amntPreAverage,'maxLocation',maxLocation,'bracketRange',bracketRange,...
            'expThreshDiffCut',expThreshDiffCut,'expThreshVoltageCut',expThreshVoltageCut,'chanInt',chanInt,'minPts',minPts,'minClustSize',minClustSize,'outlierThresh',outlierThresh);
        
    case 'average'
        [processedSig,templateArrayCellOutput] = analyFunc.template_average(templateArrayCell,rawSig,'plotIt',plotIt...,
            ,'goodVec',goodVec,'startInds',startInds,'endInds',endInds);
        
    case 'trial'
        [processedSig,templateArrayCellOutput] = analyFunc.template_trial(templateTrial,rawSig,'plotIt',plotIt...,
            ,'goodVec',goodVec,'startInds',startInds,'endInds',endInds);
end

fprintf(['-------Extracting data-------- \n \n'])

end