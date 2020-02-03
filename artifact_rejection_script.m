% David J. Caldwell - This is a script to demonstrate different stimulation artifact
% approaches to extract neural signals of interest.

% load in the data file of interest

% variables required for functions to work properly
%
% dataInt = time x channels x epochs
%
% fsData = sampling rate of the data Hz
%
% stimChans = the channels used for stimulation . These should be noted and
% excluded from further analysis
%
% plotIt = determines whether or not to plot the intermediate results of
% the functions.
%
% tEpoch = epoched time window (s)
%
% Further variables are described as they are introduced in the script
% below, as well as in the various functions called 
%
% clear the workspace
close all;clear all;clc

%%
% choose data file of interest
for dataChoice = [1]

    switch dataChoice

        case 1
            load('+data/693ffd_exampData_400ms.mat') % response timing data set
            trainDuration = [0 400]; % this is how long the stimulation train was
            xlims = [-200 1000]; % these are the x limits to visualize in plots
            chanIntList = [4 12 21 28 19 18 36 38 44 43 30 33 41 34]; % these are the channels of interest to visualize in closer detail
            minDuration = 0.5; % minimum duration of artifact in ms

        case 2
            load('+data/a1355e_examplePriming_Prime_high.mat')
            trainDuration = [0 200]; % this is how long the stimulation train was
            xlims = [-50 500]; % these are the x limits to visualize in plots
            chanIntList = [7 8 10 15 17 22 23 29 30 31 32]; % these are the channels of interest to visualize in closer detail
            minDuration = 0.5; % minimum duration of artifact in ms
            fsData = fs_data;
            tEpoch = t_epoch;
            dataInt = 4*dataInt; % needed to be multiplied by 4 from raw recording

        case 3
            load('+data/50ad9_paramSweep4.mat') % DBS data set
            fsData = fs_data;
            tEpoch = t_epoch;
            xlims = [-200 600];
            chanIntList = [5,6,7,9,10];
            trainDuration = [0 500];
            minDuration = 0.250; % minimum duration of artifact in ms
            dataInt = 4*dataInt; % needed to be multiplied by 4 from raw recording
        case 4
            load('+data/ecb43e_RHI_async_trial14.mat') % rubber hand illusion data set
            minDuration = 0.5; % minimum duration of artifact in ms
            fsData = fs_data;
            tEpoch = t_epoch;
            dataInt = 4*dataInt; % needed to be multiplied by 4 from raw recording due to acquisition parameters
            trainDuration = [0 500];
    end


    %% template subtraction
    % this is a section illustrating the template subtraction method
    % The type variable to the function call is what determines the
    % The three options
    % are 'average', 'trial', and 'dictionary'. Average is the simplest
    % approach, and on a channel by channel basis it simply averages the
    % artifact for each channel. Trial uses the stimulation pulse train within
    % each trial. 'dictionary' is the most advanced, and uses a template
    % matching algorithm with DBSCAN clustering to discover the family of
    % artifacts.
    % 'pre' defines how far back from the onset of the stimulus detect to look
    % back in time, in ms
    % post is the equivalent for after the offset of the stimulus pulse
    % fsData is the sampling frequency in Hz

    if dataChoice == 1
        type = 'dictionary';

        % if wanting to use a fixed distance, rather than dynamically detecting them
        useFixedEnd = 0;
        fixedDistance = 4; % in ms
        plotIt = 0;
        
        % parameters for detecting artifact onset and offset 
        pre = 0.8; % default time window to extend before the artifact pulse to ensure the artifact is appropriately detected (0.8 ms as default)
        post = 1; % default time window to extend before the artifact pulse to ensure the artifact is appropriately detected (1 ms as default)
        onsetThreshold = 1.5; %This value is used as absolute valued z-score threshold to determine the onset of artifacts within a train. The differentiated smoothed signal is used to determine artifact onset. This is also used in determining how many stimulation pulses are within a train, by ensuring that beginning of each artifact is within a separate artifact pulse.  
        threshVoltageCut = 75; %This is used to help determine the end of each individual artifact pulse dynamically. More specifically, this is a percentile value, against which the absolute valued, z-scored smoothed raw signal is compared to find the last value which exceeds the specified percentile voltage value. Higher values of this (i.e. closer to 100) result in a shorter duration artifact, while lower values result in a longer calculated duration of artifact. This parameter therefore should likely be set higher for more transient artifacts and lower for longer artifacts.
        threshDiffCut = 75; %This is used to help determine the end of each individual artifact pulse dynamically. More specifically, this is a percentile value, against which the absolute valued, z-scored differentiated smoothed raw signal is compared to find the last value which exceeds the specified percentile voltage value. Higher values of this (i.e. closer to 100) result in a shorter duration artifact, while lower values result in a longer calculated duration of artifact. This parameter therefore should likely be set higher for more transient artifacts and lower for longer artifacts.
        
        % these are the metrics used if the dictionary method is selected. The
        % options are 'eucl', 'cosine', 'corr', for either euclidean distance,
        % cosine similarity, or correlation for clustering and template matching.
        distanceMetricDbscan = 'eucl';
        distanceMetricSigMatch = 'corr';
        
        amntPreAverage = 3; % number of samples at the beginning of each artifact pulse to use as a baseline normalization
        normalize = 'preAverage'; % method to use for normalization of each artifact pulse. 'preAverage' uses the average across the number of samples specified above 

        % additional HDBSCAN parameters and window selection
        bracketRange = [-6:6]; %This variable sets the number of samples around the maximum voltage deflection to use for template clustering and subsequent matching. The smaller this range, the lower the dimensionality used for clustering, and the fewer points used to calculate the best matching template. This value is used to try and ensure that non-informative points are not included in the clustering and template matching. This should be set to what looks like the approximate length of the artifact's largest part.
        chanInt = 28; % channel of interest for plotting
        minPts = 2;  % Defined as k in the manuscript. This is a parameter that determines how many neighbors are used for core distance calculations for each point in the artifact window. This is a parameter that determines how many neighbors are used for core distance calculations. Increasing this parameter restricts clusters to increasingly dense areas.
        minClustSize = 3; % Defined as n in the manuscript. The minimum number of clustered artifact pulses for a cluster to be labelled as a true cluster. Increasing this number can reduce the number of clusters, and merges some clusters together that would have otherwise been labelled as individual clusters.
        outlierThresh = 0.95; % Outlier parameter for labeling artifact pulses as noise in the HDBSCAN clustering. Any artifact pulse with an outlier score greater than this will be labelled as noise. Increasing this value results in fewer points being labelled as noise 

        % optional parameters for exponential recovery, not currently used. Could be helpful for signals with large exponential recoveries
        recoverExp = 0;
        expThreshVoltageCut = 95;
        expThreshDiffCut = 95;


    elseif dataChoice == 2
        type = 'dictionary';

        useFixedEnd = 0;
        fixedDistance = 4; % in ms
        plotIt = 0;

        % parameters for detecting the onset and offset of artifacts
        pre = 0.8; % started with 1
        post = 1; % started with 0.2
        onsetThreshold = 1.5;
        threshVoltageCut = 75;
        threshDiffCut = 75;


        % these are the metrics used if the dictionary method is selected. The
        % options are 'eucl', 'cosine', 'corr', for either euclidean distance,
        % cosine similarity, or correlation for clustering ('distanceMetricDbscan') and template matching ('distanceMetricSigMatch').
        distanceMetricDbscan = 'eucl';
        distanceMetricSigMatch = 'corr';
        amntPreAverage = 3; 
        normalize = 'preAverage'; 

        bracketRange = [-6:6];
        chanInt = 15;
        minPts = 2;
        minClustSize = 3;
        outlierThresh = 0.95;

        % optional parameters for exponential recovery, not currently used. Could be helpful for signals with large exponential recoveries
        recoverExp = 0;
        expThreshVoltageCut = 95;
        expThreshDiffCut = 95;

    elseif dataChoice == 3
        type = 'dictionary';

        useFixedEnd = 0;
        fixedDistance = 4; % in ms
        plotIt = 0;

        % parameters for detecting the onset and offset of artifacts
        pre = 0.8;
        post = 1;
        onsetThreshold = 1.5;
        threshVoltageCut = 55;
        threshDiffCut = 55;

        % these are the metrics used if the dictionary method is selected. The
        % options are 'eucl', 'cosine', 'corr', for either euclidean distance,
        % cosine similarity, or correlation for clustering and template matching.

        distanceMetricDbscan = 'eucl';
        distanceMetricSigMatch = 'corr';
        amntPreAverage = 3;
        normalize = 'preAverage';

        bracketRange = [-6:6];
        chanInt = 10;
        minPts = 2;
        minClustSize = 3;
        outlierThresh = 0.95;

        % optional parameters for exponential recovery, not currently used. Could be helpful for signals with large exponential recoveries
        recoverExp = 0;
        expThreshVoltageCut = 85;
        expThreshDiffCut = 85;

    elseif dataChoice == 4
        type = 'dictionary';

        useFixedEnd = 0;
        fixedDistance = 4;
        plotIt = 0;

        pre = 0.8;
        post = 1;
        onsetThreshold = 1.5;
        threshVoltageCut = 75;
        threshDiffCut = 75;

        % these are the metrics used if the dictionary method is selected. The
        % options are 'eucl', 'cosine', 'corr', for either euclidean distance,
        % cosine similarity, or correlation for clustering and template matching.
        distanceMetricDbscan = 'eucl';
        distanceMetricSigMatch = 'corr';
        amntPreAverage = 3;
        normalize = 'preAverage';

        bracketRange = [-3:3];
        chanInt = 55;
        minPts = 5;
        minClustSize = 6;
        outlierThresh = 0.95;

        recoverExp = 0;
        expThreshVoltageCut = 95;
        expThreshDiffCut = 95;

    else
        type = 'dictionary';

        useFixedEnd = 0;
        fixedDistance = 4;
        plotIt = 0;

        pre = 0.8;
        post = 1;
        onsetThreshold = 1.5;
        threshVoltageCut = 75;
        threshDiffCut = 75;

        % these are the metrics used if the dictionary method is selected. The
        % options are 'eucl', 'cosine', 'corr', for either euclidean distance,
        % cosine similarity, or correlation for clustering and template matching.
        distanceMetricDbscan = 'eucl';
        distanceMetricSigMatch = 'corr';
        amntPreAverage = 3;
        normalize = 'preAverage';

        bracketRange = [-6:6];
        minPts = 2;
        minClustSize = 3;
        outlierThresh = 0.95;

        recoverExp = 0;
        expThreshVoltageCut = 95;
        expThreshDiffCut = 95;

    end

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
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    average = 1;
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
    end
    %
    reref = 0;
    if reref
        processedSigReref = analyFunc.rereference_CAR_median(processedSig,rerefMode,badChannels,[],[],channelsToUse);
        vizFunc.small_multiples_time_series(processedSigReref,tEpoch,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',ylims,'modePlot',modePlot,'highlightRange',trainDuration)
        fprintf(['-------Done rereferencing-------- \n'])

    end
    %     %
    %%%%%% wavelet
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

        for chanInt = chanIntList
            vizFunc.visualize_wavelet_channel_no_raw(normalizedData,tMorlet,fMorlet,processedSig,...
                tEpoch,chanInt,individual,average)
        end

        for chanInt = chanIntList
            vizFunc.visualize_wavelet_channel(normalizedData,tMorlet,fMorlet,processedSig,...
                tEpoch,dataInt,chanInt,individual,average)
        end
        %         %
        %         for chanInt = chanIntList
        %             vizFunc.visualize_wavelet_channel_small(normalizedData,tMorlet,fMorlet,processedSig,...
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
            vizFunc.visualize_raw_vs_processed(processedSig,tEpoch,dataInt,chanInt,individual,average)
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
pre = 0.8; % in ms

% this is how far to look after the algorithm detects each stimulation
% pulse onset to allow for maximal artifact rejection
post = 1; % in ms

% This is the maximal duration of time allowed for the artifact rejection
% for each pulse, and if using "useFixedEnd", it simply considers this time
% window. Otherwise, the algorithm detects the individual offset of each
% stimulation pulse.
fixedDistance = 1.2;
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
pre = 0.8;
post = 1;

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
for chanInt = chanIntList
    vizFunc.visualize_wavelet_channel_small(normalizedData,tMorlet,fMorlet,processedSig,...
        tEpoch,dataInt,chanInt,individual,average)
end
%
ylimsSpect = [5 300];
xlims = [-200 1000];
vizFunc.small_multiples_spectrogram(normalizedData,tMorlet,fMorlet,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',ylimsSpect);

%% ica example
% for the paper, it was the first data set, the chanInt = 28, the
% chanIntList = [38];

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

lnReject = false;
lnFreq = 200;
hp = false;
hpFreq = [];
lp = true;
lpFreq = [100]
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
