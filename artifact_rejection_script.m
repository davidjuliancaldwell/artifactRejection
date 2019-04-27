% David J. Caldwell - This is a script to demonstrate different stimulation artifact
% approaches to extract neural signals of interest.

% load in the data file of interest

% variables required for functions to work properly
%
% dataInt = time x channels x epochs
%
% fsData = sampling rate of the data Hz
%
% stimChans - the channels used for stimulation . These should be noted and
% excluded from further analysis
%
% plotIt - determines whether or not to plot the intermediate results of
% the functions.
%
% t_epoch - epoched time window
%%
% clear the workspace
%close all;clear all;clc
close all;clear all;clc
%%
% choose data file of interest
for dataChoice = [1]
    
    switch dataChoice
        
        case 1
            load('+data/693ffd_exampData_400ms.mat') % response timing data set
            trainDuration = [0 400]; % this is how long the stimulation train was
            xlims = [-200 1000]; % these are the x limits to visualize in plots
            chanIntList = [4 12 21 28 19 18 36 44 43 30 33 41 34]; % these are the channels of interest to visualize in closer detail
            %chanIntList = [21 28];
            minDuration = 0.5; % minimum duration of artifact in ms
            % dataInt = 4*dataInt;
            
        case 2
            load('+data/a1355e_examplePriming_Prime_high.mat')
            trainDuration = [0 200]; % this is how long the stimulation train was
            xlims = [-50 500]; % these are the x limits to visualize in plots
            chanIntList = [7 8 10 15 17 22 23 29 30 31 32]; % these are the channels of interest to visualize in closer detail
            minDuration = 0.5; % minimum duration of artifact in ms
            fsData = fs_data;
            tEpoch = t_epoch;
            %   dataInt = 4*dataInt;
            
        case 3
            load('+data/50ad9_paramSweep4.mat') % DBS data set
            fsData = fs_data;
            tEpoch = t_epoch;
            xlims = [-200 600];
            chanIntList = [5,6,7,9,10];
            trainDuration = [0 500];
            minDuration = 0.250; % minimum duration of artifact in ms
            dataInt = 4*dataInt;
        case 4
            load('+data/ecb43e_RHI_async_trial14.mat') % rubber hand illusion data set
            minDuration = 0.5; % minimum duration of artifact in ms
            fsData = fs_data;
            tEpoch = t_epoch;
            dataInt = 4*dataInt;
            trainDuration = [0 500];
    end
    
    
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
    
    if dataChoice == 1
        type = 'dictionary';
        
        useFixedEnd = 0;
        %fixedDistance = 2;
        fixedDistance = 4; % in ms
        plotIt = 0;
        
        %pre = 0.4096; % in ms
        %post = 0.4096; % in ms
        
        pre = 0.8; % started with 1
        post = 1; % started with 0.2
        % 2.8, 1, 0.5 was 3/19/2018
        
        % these are the metrics used if the dictionary method is selected. The
        % options are 'eucl', 'cosine', 'corr', for either euclidean distance,
        % cosine similarity, or correlation for clustering and template matching.
        distanceMetricDbscan = 'eucl';
        distanceMetricSigMatch = 'corr';
        amntPreAverage = 3;
        normalize = 'preAverage';
        %normalize = 'firstSamp';
        
        onsetThreshold = 1.5;
        
        recoverExp = 0;
        threshVoltageCut = 75;
        threshDiffCut = 75;
        expThreshVoltageCut = 95;
        expThreshDiffCut = 95;
        bracketRange = [-6:6];
        chanInt = 28;
        minPts = 2;
        minClustSize = 3;
        outlierThresh = 0.95;
        
        
    elseif dataChoice == 2
        type = 'dictionary';
        
        useFixedEnd = 0;
        %fixedDistance = 2;
        fixedDistance = 4; % in ms
        plotIt = 0;
        
        %pre = 0.4096; % in ms
        %post = 0.4096; % in ms
        
        pre = 0.8; % started with 1
        post = 1; % started with 0.2
        % 2.8, 1, 0.5 was 3/19/2018
        
        % these are the metrics used if the dictionary method is selected. The
        % options are 'eucl', 'cosine', 'corr', for either euclidean distance,
        % cosine similarity, or correlation for clustering and template matching.
        distanceMetricDbscan = 'eucl';
        distanceMetricSigMatch = 'corr';
        amntPreAverage = 3;
        normalize = 'preAverage';
        %normalize = 'firstSamp';
        onsetThreshold = 1.5;
        recoverExp = 0;
        threshVoltageCut = 75;
        threshDiffCut = 75;
        expThreshVoltageCut = 95;
        expThreshDiffCut = 95;
        bracketRange = [-6:6];
        chanInt = 15;
        minPts = 2;
        minClustSize = 3;
        outlierThresh = 0.95;
        
        
    elseif dataChoice == 3
        type = 'dictionary';
        
        useFixedEnd = 0;
        %fixedDistance = 2;
        fixedDistance = 4; % in ms
        plotIt = 0;
        
        %pre = 0.4096; % in ms
        %post = 0.4096; % in ms
        
        pre = 0.8; % started with 1
        post = 1; % started with 0.2
        % 2.8, 1, 0.5 was 3/19/2018
        
        % these are the metrics used if the dictionary method is selected. The
        % options are 'eucl', 'cosine', 'corr', for either euclidean distance,
        % cosine similarity, or correlation for clustering and template matching.
        
        distanceMetricDbscan = 'eucl';
        distanceMetricSigMatch = 'corr';
        amntPreAverage = 3;
        normalize = 'preAverage';
        %normalize = 'firstSamp';
        onsetThreshold = 1.5;
        recoverExp = 1;
        threshVoltageCut = 55;
        threshDiffCut = 55;
        expThreshVoltageCut = 95;
        expThreshDiffCut = 95;
        bracketRange = [-6:6];
        chanInt = 10;
        minPts = 2;
        minClustSize = 3;
        outlierThresh = 0.95;
        
    elseif dataChoice == 4
        type = 'dictionary';
        
        useFixedEnd = 0;
        %fixedDistance = 2;
        fixedDistance = 4; % in ms
        plotIt = 0;
        
        %pre = 0.4096; % in ms
        %post = 0.4096; % in ms
        
        pre = 0.8; % started with 1
        post = 1; % started with 0.2
        % 2.8, 1, 0.5 was 3/19/2018
        
        % these are the metrics used if the dictionary method is selected. The
        % options are 'eucl', 'cosine', 'corr', for either euclidean distance,
        % cosine similarity, or correlation for clustering and template matching.
        
        distanceMetricDbscan = 'eucl';
        distanceMetricSigMatch = 'corr';
        amntPreAverage = 3;
        normalize = 'preAverage';
        %normalize = 'firstSamp';
        onsetThreshold = 1.5;
        
        recoverExp = 0;
        threshVoltageCut = 75;
        threshDiffCut = 75;
        expThreshVoltageCut = 95;
        expThreshDiffCut = 95;
        bracketRange = [-3:3];
        chanInt = 55;
        minPts = 5;
        minClustSize = 6;
        outlierThresh = 0.95;
        
    else
        
        type = 'dictionary';
        
        useFixedEnd = 0;
        %fixedDistance = 2;
        fixedDistance = 4; % in ms
        plotIt = 0;
        
        %pre = 0.4096; % in ms
        %post = 0.4096; % in ms
        
        pre = 0.8; % started with 1
        post = 1; % started with 0.2
        % 2.8, 1, 0.5 was 3/19/2018
        
        % these are the metrics used if the dictionary method is selected. The
        % options are 'eucl', 'cosine', 'corr', for either euclidean distance,
        % cosine similarity, or correlation for clustering and template matching.
        distanceMetricDbscan = 'eucl';
        distanceMetricSigMatch = 'corr';
        amntPreAverage = 3;
        normalize = 'preAverage';
        %normalize = 'firstSamp';
        
        onsetThreshold = 1.5;
        
        recoverExp = 0;
        threshVoltageCut = 75;
        threshDiffCut = 75;
        expThreshVoltageCut = 95;
        expThreshDiffCut = 95;
        bracketRange = [-6:6];
        minPts = 2;
        minClustSize = 3;
        outlierThresh = 0.95;
        
        
    end
    
    
%%          
    % set parameters for fit function
    fRange = [59 61]; 
    smoothSpan = 51; 
    tic;
    [subtractedSig,phase,f,r,fitline] = analyFunc.sine_fit(dataInt,tEpoch,smoothSpan,fRange,fsData,plotIt);
toc;
        fRange = [119 121]; 
        [subtractedSig2,phase,f,r,fitline] = analyFunc.sine_fit(subtractedSig,tEpoch,smoothSpan,fRange,fsData,plotIt);

              fRange = [179 181]; 
        [subtractedSig3,phase,f,r,fitline] = analyFunc.sine_fit(subtractedSig2,tEpoch,smoothSpan,fRange,fsData,plotIt);

%%    
    [processedSig,templateDictCell,templateTrial,startInds,endInds] = analyFunc.template_subtract(dataInt,'type',type,...
        'fs',fsData,'plotIt',plotIt,'pre',pre,'post',post,'stimChans',stimChans,...
        'useFixedEnd',useFixedEnd,'fixedDistance',fixedDistance,...,
        'distanceMetricDbscan',distanceMetricDbscan,'distanceMetricSigMatch',distanceMetricSigMatch,...
        'recoverExp',recoverExp,'normalize',normalize,'amntPreAverage',amntPreAverage,...
        'minDuration',minDuration,'bracketRange',bracketRange,'threshVoltageCut',threshVoltageCut,...
        'threshDiffCut',threshDiffCut,'expThreshVoltageCut',expThreshVoltageCut,...
        'expThreshDiffCut',expThreshDiffCut,'onsetThreshold',onsetThreshold,'chanInt',chanInt,...
        'minPts',minPts,'minClustSize',minClustSize,'outlierThresh',outlierThresh);
    
    % visualization
    % of note - more visualizations are created here, including what the
    % templates look like on each channel, and what the discovered templates are
    xlims = [-100 500];
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     vizFunc.multiple_visualizations(processedSig,dataInt,'fs',fsData,'type',type,'tEpoch',...
    %         tEpoch,'xlims',xlims,'trainDuration',trainDuration,'stimChans',stimChans,...,
    %         'chanIntList',chanIntList,'templateTrial',templateTrial,'templateDictCell',templateDictCell,'modePlot','confInt')
    %     %
    average = 1;
    %chanIntList = 3;
    %  trainDuration = [0 400];
    modePlot = 'avg';
    xlims = [-200 1000];
    ylims = [-0.6 0.6];
    vizFunc.small_multiples_time_series(processedSig,tEpoch,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',ylims,'modePlot',modePlot,'highlightRange',trainDuration)
    
    %% additional processing
    
    rerefMode = 'selectedChannelsMean';
    switch dataChoice
        case 1
            badChannels = [stimChans [53:64]];
            channelsToUse = [1:8 41:48];
            channelsToUse = [49 42 35 14 7 8];
        case 2
            badChannels = stimChans;
            channelsToUse = [22 23 30 31 38 39 46 47];
        case 6
            badChannels = stimChans;
        case 7
            badChannels = stimChans;
            channels
    end
    %
    reref = 0;
    if reref
        processedSigReref = analyFunc.rereference_CAR_median(processedSig,rerefMode,badChannels,[],[],channelsToUse);
        vizFunc.small_multiples_time_series(processedSigReref,tEpoch,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',ylims,'modePlot',modePlot,'highlightRange',trainDuration)
        fprintf(['-------Done rereferencing-------- \n'])
        
    end
    %%
    %%%%%%% wavelet
    fprintf(['-------Beginning wavelet analysis-------- \n'])
    
    timeRes = 0.01; % 10 ms bins
    
    [powerout,fMorlet,tMorlet,~] = analyFunc.waveletWrapper(processedSig,fsData,timeRes,stimChans);
    %
    fprintf(['-------Ending wavelet analysis-------- \n'])
    
    % additional parameters
    postStim = 2000;
    sampsPostStim = round(postStim/1e3*fsData);
    
    preStim = 1000;
    sampsPreStim = round(preStim/1e3*fsData);
    
    tMorlet = linspace(-preStim,postStim,length(tMorlet))/1e3;
    % normalize data
    dataRef = powerout(:,tMorlet<0.05 & tMorlet>-0.8,:,:);
    %
    [normalizedData] = analyFunc.normalize_spectrogram_wavelet(dataRef,powerout);
    individual = 0;
    average = 1;
    %%
    if dataChoice == 1 || dataChoice == 2
        
        %         for chanInt = chanIntList
        %             vizFunc.visualize_wavelet_channel_no_raw(normalizedData,tMorlet,fMorlet,processedSigReref,...
        %                 tEpoch,chanInt,individual,average)
        %         end
        %
        for chanInt = chanIntList
            vizFunc.visualize_wavelet_channel(normalizedData,tMorlet,fMorlet,processedSig,...
                tEpoch,dataInt,chanInt,individual,average)
        end
        %         %
        %         for chanInt = chanIntList
        %             vizFunc.visualize_wavelet_channel_small(normalizedData,tMorlet,fMorlet,processedSigReref,...
        %                 tEpoch,dataInt,chanInt,individual,average)
        %       end
        %
        ylimsSpect = [5 300];
        xlims = [-200 1000];
        vizFunc.small_multiples_spectrogram(normalizedData,tMorlet,fMorlet,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',ylimsSpect);
    end
    %%
    if dataChoice == 3
        individual = 0;
        average = 1;
        for chanInt = chanIntList
            vizFunc.visualize_raw_vs_processed(processedSigReref,tEpoch,dataInt,chanInt,individual,average)
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
return
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
    'stimchans',stimChans,'useFixedEnd',useFixedEnd,'fixedDistance',fixedDistance,'pre',pre,'post',post,'onsetThreshold',onsetThreshold);

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
pre = 0.8; % started with 1
post = 1; % started with 0.2
% 2.8, 1, 0.5 was 3/19/2018

fixedDistance = 2; % in ms

[processedSig,startInds,endInds] = analyFunc.interpolate_artifact(dataInt,'fs',fsData,...
    'plotIt',plotIt,'type',type,'stimchans',stimChans,'useFixedEnd',useFixedEnd,...
    'fixedDistance',fixedDistance,'pre',pre,'post',post,'onsetThreshold',onsetThreshold,...
    'threshVoltageCut',threshVoltageCut,'threshDiffCut',threshDiffCut);

%visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vizFunc.multiple_visualizations(processedSig,dataInt,'fs',fsData,'type',type,'tEpoch',...
    tEpoch,'xlims',xlims,'trainDuration',trainDuration,'stimChans',stimChans,...,
    'chanIntList',chanIntList,'modePlot','confInt')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
average = 1;
%chanIntList = 3;
trainDuration = [0 400];
modePlot = 'avg';
xlims = [-200 1000];
ylims = [-0.6 0.6];
vizFunc.small_multiples_time_series(processedSig,tEpoch,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',ylims,'modePlot',modePlot,'highlightRange',trainDuration)

rerefMode = 'mean';
switch dataChoice
    case 1
        badChannels = [stimChans [53:64]];
    case 2
        badChannels = stimChans;
        channelsToUse = [22 23 30 31 38 39 46 47];
    case 6
        badChannels = stimChans;
    case 7
        badChannels = stimChans;
end
reref = 0;
if reref
    processedSig = analyFunc.rereference_CAR_median(processedSig,rerefMode,badChannels,[],[],channelsToUse);
    processedSigReref = analyFunc.rereference_CAR_median(processedSig,rerefMode,badChannels,[],[],channelsToUse);
end

%
%%%%%%% wavelet
timeRes = 0.01; % 25 ms bins

% [powerout,fMorlet,tMorlet] = wavelet_wrapper(processedSig,fsData,stimChans);
[powerout,fMorlet,tMorlet,~] = analyFunc.waveletWrapper(processedSig,fsData,timeRes,stimChans);
%
% additional parameters
postStim = 2000;
sampsPostStim = round(postStim/1e3*fsData);

preStim = 1000;
sampsPreStim = round(preStim/1e3*fsData);

tMorlet = linspace(-preStim,postStim,length(tMorlet))/1e3;
% normalize data
dataRef = powerout(:,tMorlet<0.05 & tMorlet>-0.8,:,:);
%
[normalizedData] = analyFunc.normalize_spectrogram_wavelet(dataRef,powerout);
individual = 0;
average = 1;
% chanIntLIst = 42;
%     %%
%     % chanIntList = chanInt;
%     for chanInt = chanIntList
%         vizFunc.visualize_wavelet_channel_no_raw(normalizedData,tMorlet,fMorlet,processedSig,...
%             tEpoch,chanInt,individual,average)
%     end
%
%      for chanInt = chanIntList
%         vizFunc.visualize_wavelet_channel(normalizedData,tMorlet,fMorlet,processedSig,...
%             tEpoch,dataInt,chanInt,individual,average)
%      end
%     %%
%       for chanInt = chanIntList
%         vizFunc.visualize_wavelet_channel_small(normalizedData,tMorlet,fMorlet,processedSig,...
%             tEpoch,dataInt,chanInt,individual,average)
%     end
%
ylimsSpect = [5 300];
xlims = [-200 1000];
vizFunc.small_multiples_spectrogram(normalizedData,tMorlet,fMorlet,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',ylimsSpect);

%% ica example

chanInt = 28;
scaleFactor = 600;
numComponentsSearch = 20;
plotIt = false;
meanSub = 0;
orderPoly = 3;
[processedSig,subtractedSigCell,reconArtifactMatrix,reconArtifact,t] = analyFunc.ica_artifact_remove_train(tEpoch,dataInt,stimChans,fsData,scaleFactor,numComponentsSearch,plotIt,chanInt,meanSub,orderPoly);
%
%%%%%%% wavelet
timeRes = 0.01; % 25 ms bins

% [powerout,fMorlet,tMorlet] = wavelet_wrapper(processedSig,fsData,stimChans);
[powerout,fMorlet,tMorlet,~] = analyFunc.waveletWrapper(processedSig,fsData,timeRes,stimChans);
%
% additional parameters
postStim = 2000;
sampsPostStim = round(postStim/1e3*fsData);

preStim = 1000;
sampsPreStim = round(preStim/1e3*fsData);

tMorlet = linspace(-preStim,postStim,length(tMorlet))/1e3;
% normalize data
dataRef = powerout(:,tMorlet<0.05 & tMorlet>-0.8,:,:);
%
[normalizedData] = analyFunc.normalize_spectrogram_wavelet(dataRef,powerout);
individual = 0;
average = 1;
%%
chanIntList = [38];
% chanIntList = chanInt;
for chanInt = chanIntList
    vizFunc.visualize_wavelet_channel_no_raw(normalizedData,tMorlet,fMorlet,processedSig,...
        tEpoch,chanInt,individual,average)
end
%
for chanInt = chanIntList
    vizFunc.visualize_wavelet_channel(normalizedData,tMorlet,fMorlet,processedSig,...
        tEpoch,dataInt,chanInt,individual,average)
end
%
for chanInt = chanIntList
    vizFunc.visualize_wavelet_channel_small(normalizedData,tMorlet,fMorlet,processedSig,...
        tEpoch,dataInt,chanInt,individual,average)
end

%%
% visualization
% of note - more visualizations are created here, including what the
% templates look like on each channel, and what the discovered templates are
xlims = [-100 500];
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     vizFunc.multiple_visualizations(processedSig,dataInt,'fs',fsData,'type',type,'tEpoch',...
%         tEpoch,'xlims',xlims,'trainDuration',trainDuration,'stimChans',stimChans,...,
%         'chanIntList',chanIntList,'templateTrial',templateTrial,'templateDictCell',templateDictCell,'modePlot','confInt')
%     %
average = 1;
%chanIntList = 3;
%  trainDuration = [0 400];
modePlot = 'avg';
xlims = [-200 1000];
ylims = [-0.6 0.6];
vizFunc.small_multiples_time_series(processedSig,tEpoch,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',ylims,'modePlot',modePlot,'highlightRange',trainDuration)

%% low pass filter

% lnReject = false;
% lnFreq = [];
% hp = false;
% hpFreq = [];
% lp = true;
% lpFreq = 100;
% filterOrder = 4;
% causality = 'acausal';

lnReject = false;
lnFreq = 200;
hp = false;
hpFreq = [];
lp = true;
lpFreq = [100];p
filterOrder = 4;
causality = 'acausal';
for index = 1:size(dataInt,3)
    processedSig(:,:,index) = ecogFilter(squeeze(dataInt(:,:,index)), lnReject, lnFreq, hp, hpFreq, lp, lpFreq, fsData, filterOrder, causality);
end

%
%%%%%%% wavelet
timeRes = 0.01; % 25 ms bins

% [powerout,fMorlet,tMorlet] = wavelet_wrapper(processedSig,fsData,stimChans);
[powerout,fMorlet,tMorlet,~] = analyFunc.waveletWrapper(processedSig,fsData,timeRes,stimChans);
%
% additional parameters
postStim = 2000;
sampsPostStim = round(postStim/1e3*fsData);

preStim = 1000;
sampsPreStim = round(preStim/1e3*fsData);

tMorlet = linspace(-preStim,postStim,length(tMorlet))/1e3;
% normalize data
dataRef = powerout(:,tMorlet<0.05 & tMorlet>-0.8,:,:);
%
[normalizedData] = analyFunc.normalize_spectrogram_wavelet(dataRef,powerout);
individual = 0;
average = 1;
%
chanIntList = 28;
% chanIntList = chanInt;
for chanInt = chanIntList
    vizFunc.visualize_wavelet_channel_no_raw(normalizedData,tMorlet,fMorlet,processedSig,...
        tEpoch,chanInt,individual,average)
end
%
for chanInt = chanIntList
    vizFunc.visualize_wavelet_channel(normalizedData,tMorlet,fMorlet,processedSig,...
        tEpoch,dataInt,chanInt,individual,average)
end
%
for chanInt = chanIntList
    vizFunc.visualize_wavelet_channel_small(normalizedData,tMorlet,fMorlet,processedSig,...
        tEpoch,dataInt,chanInt,individual,average)
end
