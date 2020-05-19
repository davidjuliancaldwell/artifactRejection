%% make haptic touch plots for manuscript
close all; clearvars ; clc

sid = '693ffd';
reref = 0;
rerefMode = 'mean';
primedBlock = 0;
type = 'linear';

load('+data/693ffdpooledData_changePts_noDelay.mat') % response timing data set

if ~exist('fsData','var')
    fsData = fs_data;
end
%% combine the pooled data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('tEpoch','var')
    tEpoch = t_epoch;
end

[buttonLocsSamps,buttonLocs,data,tEpoch,uniqueCond] = analyFunc.combine_pooled_data(sid,epochedCortEco_cell,tEpoch);

% additional parameters
postStim = 2000;
sampsPostStim = round(postStim/1e3*fsData);

preStim = 1000;
sampsPreStim = round(preStim/1e3*fsData);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

chanInt = 28;
stimChans = [20 29];
chanIntList = [21 28 19 36 44 43 30];

%%
% condIntAns for RT task
% -1 = tactor
%  0 = null
%  1 = off-target
%  2-5 = conditions 1->4

% condIntAns for priming
% 0 - unprimed
% 1 - primed

% for 3ada8b, block 2
%
% uniqueCond =
%
%      0
%      5
%      6 - this is the 2 pulses in isolation

condInt = 1;

condIntAns = uniqueCond(condInt);
dataInt = data{condInt};

if ~strcmp(sid,'3ada8b')
    dataInt = 4*dataInt;
end
response = buttonLocs{condInt};

% get additional fake channels on 693ffd so it plots ok
if strcmp(sid,'693ffd')
    dataInt = cat(2,dataInt,zeros(size(dataInt,1),64-size(dataInt,2),size(dataInt,3)));
    stimChans = [stimChans, [53:64]];
end

buttonLocsInt = buttonLocs{condInt};
processedSig = dataInt;

%stimTime = 1e3*tactorLocsVec; %
stimTime = zeros(size(processedSig,3),1); % it is centered around zero now
t = tEpoch;

% this selects trials with all response times
responseBool = logical(ones(size(response)));

response = response(responseBool);
% only take ones where they responded within time bins
if ~strcmp(sid,'3ada8b') & (primedBlock == 2);
    processedSig = processedSig(:,:,responseBool);
end

badChannels = [stimChans [53:64]];

if reref
    processedSig = analyFunc.rereference_CAR_median(processedSig,rerefMode,badChannels);
end

processedSigReref = analyFunc.rereference_CAR_median(processedSig,rerefMode,badChannels);

% small multiples

individual = 0;
average = 1;
%chanIntList = 3;
trainDuration = [];
modePlot = 'avg';
xlims = [-200 1000];
ylims = [-300 300];
vizFunc.small_multiples_time_series(processedSig,tEpoch,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',ylims,'modePlot',modePlot,'highlightRange',trainDuration)
vizFunc.small_multiples_time_series(processedSigReref,tEpoch,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',ylims,'modePlot',modePlot,'highlightRange',trainDuration)

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trainDuration = [];
vizFunc.multiple_visualizations(processedSig,dataInt,'fs',fsData,'type',type,'tEpoch',...
    tEpoch,'xlims',xlims,'trainDuration',trainDuration,'stimChans',stimChans,...,
    'chanIntList',chanIntList,'modePlot','confInt')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% wavelet
timeRes = 0.01; % 10 ms bins

[powerout,fMorlet,tMorlet,~] = analyFunc.waveletWrapper(processedSig,fsData,timeRes,stimChans);
%
tMorlet = linspace(-preStim,postStim,length(tMorlet))/1e3;
% normalize data
dataRef = powerout(:,tMorlet<0.05 & tMorlet>-0.8,:,:);
%
[normalizedData] = analyFunc.normalize_spectrogram_wavelet(dataRef,powerout);

[normalizedData_avg] = analyFunc.normalize_spectrogram_wavelet_avg(dataRef,powerout);

[normalizedData_baseline_fusion] = analyFunc.normalize_spectrogram_wavelet_fused_baseline(dataRef,powerout);

normalizedData_baseline_on_avg = analyFunc.normalize_spectrogram_wavelet_on_avg_signal(dataRef,powerout);

individual = 0;
average = 1;
% chanIntLIst = 42;

%%
[poweroutReref,fMorlet,tMorlet,~] = analyFunc.waveletWrapper(processedSigReref,fsData,timeRes,stimChans);
tMorlet = linspace(-preStim,postStim,length(tMorlet))/1e3;

dataRefReref = poweroutReref(:,tMorlet<0.05 & tMorlet>-0.8,:,:);

[normalizedDataReref] = analyFunc.normalize_spectrogram_wavelet(dataRefReref,poweroutReref);

%%
% chanIntList = chanInt;
for chanInt = chanIntList
    %visualize_wavelet_channel(normalizedData,tMorlet,fMorlet,processedSig,...
    %  tEpoch,dataInt,chanInt,stimTime,response,individual,average)
    
    vizFunc.visualize_wavelet_channel(normalizedData,tMorlet,fMorlet,processedSig,...
        tEpoch,dataInt,chanInt,individual,average,xlims)
    
    vizFunc.visualize_wavelet_channel(normalizedData_avg,tMorlet,fMorlet,processedSig,...
        tEpoch,dataInt,chanInt,individual,average,xlims)
    
    vizFunc.visualize_wavelet_channel(normalizedData_baseline_fusion,tMorlet,fMorlet,processedSig,...
        tEpoch,dataInt,chanInt,individual,average,xlims)
    
        vizFunc.visualize_wavelet_channel(normalizedData_baseline_on_avg,tMorlet,fMorlet,processedSig,...
        tEpoch,dataInt,chanInt,individual,average,xlims)
    
   % vizFunc.visualize_wavelet_channel(normalizedDataReref,tMorlet,fMorlet,processedSigReref,...
   %     tEpoch,dataInt,chanInt,individual,average,xlims)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

HGPowerWavelet = squeeze(mean(squeeze(normalizedDataReref(fMorlet < 150 & fMorlet > 70,:,:,:)),1));

%%
vizFunc.small_multiples_spectrogram(normalizedDataReref,tMorlet,fMorlet,'type1',stimChans,'type2',0,'xlims',xlims);
%% hilb amp HG
processedSigHG = zeros(size(processedSig));
for trial = 1:size(processedSig,3)
    [amp] = log(hilbAmp(squeeze(processedSig(:,:,trial)), [70 150], fsData).^2);
    processedSigHG(:,:,trial) = amp;
end
%%
%chanInt = 1;
figure
subplot(2,1,1)
plot(1e3*tEpoch,squeeze(mean(squeeze(processedSigHG(:,chanInt,:)),2)))
xlabel('time (ms)')
ylabel('power (log(HG amplitude squared)')
xlim([-50 500])
vline(0)
set(gca,'fontsize',14)

title(['hilbert HG amplitude - channel ' num2str(chanInt)])
subplot(2,1,2)
plot(1e3*tMorlet,mean(squeeze(HGPowerWavelet(:,chanInt,:)),2))
xlim([-50 500])
vline(0)
xlabel('time (ms)')
ylabel('power normalized to baseline')
title(['average wavelet amplitude - channel ' num2str(chanInt)])
set(gca,'fontsize',14)
%%
figure
trials = 1:size(HGPowerWavelet,3);
time = tMorlet;
tLow = -0.2;
tHigh = 0.5;
imagesc(1e3*tMorlet(tMorlet>tLow & tMorlet < tHigh),trials,squeeze(HGPowerWavelet((tMorlet>tLow & tMorlet < tHigh),chanInt,:))')
colormap(flipud(bone))
axis('normal')
ylabel('trial')
xlabel('time (ms)')
colorbar()
title('average wavelet HG amplitude')
set(gca,'fontsize',14)
%%
%vizFunc.small_multiples_time_series(processedSigHG,tEpoch,'type1',stimChans,'type2',0,'xlims',xlims,'ylims',[-40 -20],'modePlot','avg','highlightRange',trainDuration)


