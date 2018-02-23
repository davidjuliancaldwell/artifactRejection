% David Caldwell - This is a script to demonstrate different stimulation artifact
% approaches to extract neural signals of interest.

% clear the workspace
close all;clear all;clc

% load in the data file of interest
dataChoice = 1;

switch dataChoice
        case 1
        load('693ffd_exampData_800ms.mat') % response timing data set
        train_duration = [0 800]; % this is how long the stimulation train was
        xlims = [-200 1000]; % these are the x limits to visualize in plots
        chanIntList = [12 21 28 19 18 36 44 43 30 33 41 34]; % these are the channels of interest to visualize in closer detail
    case 2
        load('2fd831_exampData_400ms.mat') % response timing data set
        train_duration = [0 400]; % this is how long the stimulation train was
        xlims = [-200 1000]; % these are the x limits to visualize in plots
        chanIntList = [2 10 51 42]; % these are the channels of interest to visualize in closer detail
    case 3
        load('50ad9_paramSweep4.mat') % DBS data set
    case 4
        load('ecb43e_RHI_async_trial14.mat') % rubber hand illusion data set
    case 5
        load('3f2113_stim_12_52.mat') % stimulation spacing data set
        train_duration = [0 5];
        xlims = [-50 200];
        
end

% variables required for functions to work properly
%
% dataInt = time x channels x trials
%
% fs_data = sampling rate of the data Hz
%
% stimChans - the channels used for stimulation . These should be noted and
% exluded from further analysis
%
% plotIt - determines whether or not to plot the intermediate results of
% the functions.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% linear interpolation with simple linear interpolation scheme

% process the signal

% the type variable here determines whether to use a linear interpolation
% scheme or a polynomial spline interpolation scheme
type = 'linear';

% this determines whether or not to march a set amount of time after
% stimulation onset, or to detect the end of each pulse
useFixedEnd = 1; 

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
fixed_distance = 2.2; % in ms

% perform the processing
processedSig = interpolate_artifact(dataInt,'fs',fs_data,'plotIt',0,'type',type,...,
    'stimchans',stimChans,'useFixedEnd',useFixedEnd,'fixed_distance',fixed_distance,'pre',pre,'post',post);

% visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% multiple visualizations are created with this call, including time domain
% and spectral analyses

multiple_visualizations(processedSig,dataInt,'fs_data',fs_data,'type',type,'t_epoch',...
    t_epoch,'xlims',xlims,'train_duration',train_duration,'stimChans',stimChans,...,
    'chanIntList',chanIntList)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% linear interpolation with polynomial piecewise interpolation

% using piecewise polynomial interpolation here 
type = 'pchip';
useFixedEnd = 1;
%pre = 0.4096; % in ms
pre = 0.6;
post = 0.4096; % in ms
fixed_distance = 2.2; % in ms


processedSig = interpolate_artifact(dataInt,'fs',fs_data,'plotIt',1,'type',type,...,
    'stimchans',stimChans,'useFixedEnd',useFixedEnd,'fixed_distance',fixed_distance,'pre',pre,'post',post);

% visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

multiple_visualizations(processedSig,dataInt,'fs_data',fs_data,'type',type,'t_epoch',...
    t_epoch,'xlims',xlims,'train_duration',train_duration,'stimChans',stimChans,...,
    'chanIntList',chanIntList)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% template subtraction
% this is a section illustrating the template subtraction method
% The type variable to the function call is what determines the
% method used

% pre defines how far back from the onset of the stimulus detect to look
% back in time, in ms
% post is the equivalent for after the offset of the stimulus pulse
% fs_data is the sampling frequency in Hz

% this determines the type of template approach to use. The three options
% are 'average', 'trial', and 'dictionary'. Average is the simplest
% approach, and on a channel by channel basis it simply averages the
% artifact for each channel. Trial uses the stimulation pulse train within
% each trial. 'Dictionary' is the most advanced, and uses a template
% matching algorithm with DBSCAN clustering to discover the family of
% artifacts. 

type = 'trial';
useFixedEnd = 1;
fixed_distance = 2.8; % in ms

%pre = 0.4096; % in ms
%post = 0.4096; % in ms

pre = 1; % started with 0.7, then 0.9, adjusted the window, move it to 0.4
post = 0.2; % started with 0.7, then 1.1, then 1.5, then 1.8

% these are the metrics used if the dictionary method is selected. The
% options are 'eucl', 'cosine', 'corr', for either euclidean distance,
% cosine similarity, or correlation for clustering and template matching. 

distance_metric_dbscan = 'eucl';
distance_metric_sigMatch = 'eucl';
[processedSig,templateDict_cell,template,start_inds,end_inds] = templateSubtract(dataInt,'type',type,...
    'fs',fs_data,'plotIt',0,'pre',pre,'post',post,'stimChans',stimChans,'useFixedEnd',useFixedEnd,'fixed_distance',fixed_distance,...,
    'distance_metric_dbscan',distance_metric_dbscan,'distance_metric_sigMatch',distance_metric_sigMatch);

% visualization
% of note - more visualizations are created here, including what the
% templates look like on each channel, and what the discovered templates are
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
multiple_visualizations(processedSig,dataInt,'fs_data',fs_data,'type',type,'t_epoch',...
    t_epoch,'xlims',xlims,'train_duration',train_duration,'stimChans',stimChans,...,
    'chanIntList',chanIntList,'template',template,'templateDict_cell',templateDict_cell)

