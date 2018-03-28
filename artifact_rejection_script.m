% David Caldwell - This is a script to demonstrate different stimulation artifact
% approaches to extract neural signals of interest.

% clear the workspace
close all;clear all;clc

% load in the data file of interest

% variables required for functions to work properly
%
% dataInt = time x channels x trials
%
% fsData = sampling rate of the data Hz
%
% stimChans - the channels used for stimulation . These should be noted and
% exluded from further analysis
%
% plotIt - determines whether or not to plot the intermediate results of
% the functions.
%
% t_epoch - epoched time window


dataChoice = 6;

switch dataChoice
    
    case 1
        load('+data/693ffd_exampData_800ms.mat') % response timing data set
        trainDuration = [0 800]; % this is how long the stimulation train was
        xlims = [-200 1000]; % these are the x limits to visualize in plots
        chanIntList = [12 21 28 19 18 36 44 43 30 33 41 34]; % these are the channels of interest to visualize in closer detail
        minDuration = 0.5; % minimum duration of artifact in ms
    case 2
        load('+data/2fd831_exampData_400ms.mat') % response timing data set
        trainDuration = [0 400]; % this is how long the stimulation train was
        xlims = [-200 1000]; % these are the x limits to visualize in plots
        chanIntList = [2 10 51 42]; % these are the channels of interest to visualize in closer detail
                minDuration = 0.5; % minimum duration of artifact in ms

    case 3
        load('+data/3f2113_stim_12_52.mat') % stimulation spacing data set
        trainDuration = [0 5];
        xlims = [-10 100];
        chanIntList = [13 53 51 42 60 61 ]; % these are the channels of interest to visualize in closer detail
                minDuration = 0.5; % minimum duration of artifact in ms

    case 4
        load('+data/50ad9_paramSweep4.mat') % DBS data set
        xlims = [-200 600];
        chanIntList = [1:10];
        trainDuration = [0 500];
                minDuration = 0.250; % minimum duration of artifact in ms

    case 5
        load('+data/ecb43e_RHI_async_trial14.mat') % rubber hand illusion data set
    case 6
        load('+data/a1355e_examplePriming_Prime_high.mat')
        trainDuration = [0 200]; % this is how long the stimulation train was
        xlims = [-200 1000]; % these are the x limits to visualize in plots
        chanIntList = [7 8 10 15 17 22 23 29 30 31 32]; % these are the channels of interest to visualize in closer detail
                minDuration = 0.5; % minimum duration of artifact in ms

    case 7
        load('+data/a1355e_examplePriming_noPrime_high.mat')
        trainDuration = [0 200]; % this is how long the stimulation train was
        xlims = [-200 1000]; % these are the x limits to visualize in plots
        chanIntList = [7 8 10 15 17 22 23 29 30 31 32]; % these are the channels of interest to visualize in closer detail
                minDuration = 0.5; % minimum duration of artifact in ms

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% linear interpolation with simple linear interpolation scheme

% process the signal

% the type variable here determines whether to use a linear interpolation
% scheme or a polynomial spline interpolation scheme
type = 'linear';

% this determines whether or not to march a set amount of time after
% stimulation onset, or to detect the end of each pulse dynamically
useFixedEnd = 0;

% this is how far to look before the algorithm detects each stimulation
% pulse onset to allow for maximal artifact rejection
pre = 0.6; % in ms


% this is how far to look after the algorithm detects each stimulation
% pulse onset to allow for maximal artifact rejection
post = 0.4096; % in ms

% This is the maximal duration of time allowed for the artifact rejection
% for each pulse, and if using "useFixedEnd", it simply considers this time
% window. Otherwise, the algorithm detects the individual offset of each
% stimulation pulse.
fixedDistance = 1.2; % in ms % 2.2 for the first 2 cases, 4 for the 3rd,

% perform the processing
[processedSig,startInds,endInds] = analyFunc.interpolate_artifact(dataInt,'fs',fsData,'plotIt',0,'type',type,...,
    'stimchans',stimChans,'useFixedEnd',useFixedEnd,'fixedDistance',fixedDistance,'pre',pre,'post',post);

% visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% multiple visualizations are created with this call, including time domain
% and spectral analyses

vizFunc.multiple_visualizations(processedSig,dataInt,'fs',fsData,'type',type,'tEpoch',...
    tEpoch,'xlims',xlims,'trainDuration',trainDuration,'stimChans',stimChans,...,
    'chanIntList',chanIntList,'modePlot','confInt')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% linear interpolation with polynomial piecewise interpolation

% using piecewise polynomial interpolation here
type = 'pchip';
useFixedEnd = 0;
%pre = 0.4096; % in ms
pre = 0.6;
post = 0.4096; % in ms1
fixedDistance = 2; % in ms

[processedSig,startInds,endInds] = analyFunc.interpolate_artifact(dataInt,'fs',fsData,'plotIt',0,'type',type,...,
    'stimchans',stimChans,'useFixedEnd',useFixedEnd,'fixedDistance',fixedDistance,'pre',pre,'post',post);

% visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vizFunc.multiple_visualizations(processedSig,dataInt,'fs',fsData,'type',type,'tEpoch',...
    tEpoch,'xlims',xlims,'trainDuration',trainDuration,'stimChans',stimChans,...,
    'chanIntList',chanIntList,'modePlot','confInt')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% template subtraction
% this is a section illustrating the template subtraction method
% The type variable to the function call is what determines the
% method used

% pre defines how far back from the onset of the stimulus detect to look
% back in time, in ms
% post is the equivalent for after the offset of the stimulus pulse
% fsData is the sampling frequency in Hz

% this determines the type of template approach to use. The three options
% are 'average', 'trial', and 'dictionary'. Average is the simplest
% approach, and on a channel by channel basis it simply averages the
% artifact for each channel. Trial uses the stimulation pulse train within
% each trial. 'dictionary' is the most advanced, and uses a template
% matching algorithm with DBSCAN clustering to discover the family of
% artifacts.


if dataChoice == 6 || dataChoice == 7
    type = 'dictionary';

    useFixedEnd = 0;
    %fixedDistance = 2;
    fixedDistance = 4; % in ms
    plotIt = 0;
    
    %pre = 0.4096; % in ms
    %post = 0.4096; % in ms
    
    pre = 0.8; % started with 1
    post = 0.5; % started with 0.2
    % 2.8, 1, 0.5 was 3/19/2018
    
    % these are the metrics used if the dictionary method is selected. The
    % options are 'eucl', 'cosine', 'corr', for either euclidean distance,
    % cosine similarity, or correlation for clustering and template matching.
    
    distanceMetricDbscan = 'cosine';
    distanceMetricSigMatch = 'eucl';
    amntPreAverage = 3;
    normalize = 'preAverage';
    %normalize = 'firstSamp';
    
    recoverExp = 0;
    
elseif dataChoice == 4
    type = 'dictionary';

    useFixedEnd = 0;
    %fixedDistance = 2;
    fixedDistance = 4; % in ms
    plotIt = 0;
   
    %pre = 0.4096; % in ms
    %post = 0.4096; % in ms
    
    pre = 0.4; % started with 1
    post = 0.2; % started with 0.2
    % 2.8, 1, 0.5 was 3/19/2018
    
    % these are the metrics used if the dictionary method is selected. The
    % options are 'eucl', 'cosine', 'corr', for either euclidean distance,
    % cosine similarity, or correlation for clustering and template matching.
    
    distanceMetricDbscan = 'cosine';
    
    distanceMetricSigMatch = 'eucl';
    amntPreAverage = 3;
    normalize = 'preAverage';
    %normalize = 'firstSamp';
    
    recoverExp = 0;    
   
    
end

[processedSig,templateDictCell,templateTrial,startInds,endInds] = analyFunc.template_subtract(dataInt,'type',type,...
    'fs',fsData,'plotIt',plotIt,'pre',pre,'post',post,'stimChans',stimChans,'useFixedEnd',useFixedEnd,'fixedDistance',fixedDistance,...,
    'distanceMetricDbscan',distanceMetricDbscan,'distanceMetricSigMatch',distanceMetricSigMatch,...
    'recoverExp',recoverExp,'normalize',normalize,'amntPreAverage',amntPreAverage,'minDuration',minDuration);
%
% visualization
% of note - more visualizations are created here, including what the
% templates look like on each channel, and what the discovered templates are
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vizFunc.multiple_visualizations(processedSig,dataInt,'fs',fsData,'type',type,'tEpoch',...
    tEpoch,'xlims',xlims,'trainDuration',trainDuration,'stimChans',stimChans,...,
    'chanIntList',chanIntList,'templateTrial',templateTrial,'templateDictCell',templateDictCell,'modePlot','confInt')

%%
[processedSig_v2,templateDictCell,template] = analyFunc.template_subtract_iterative(processedSig,...,
    'fs',fsData,'plotIt',0,'pre',pre,'post',post,'stimChans',stimChans,'startInds',startInds,'endInds',endInds);
