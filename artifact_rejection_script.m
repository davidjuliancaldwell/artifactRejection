% DJC - This is a script to demonstrate different stimulation artifact
% approaches to extract neural signals of interest.

% clear the workspace
close all;clear all;clc

% load in the data file of interest
dataChoice = 1;

switch dataChoice
    case 1
        load('2fd831_exampData_400ms.mat') % response timing data set
    case 2
        load('693ffd_exampData_800ms.mat') % response timing data set
    case 3
        load('50ad9_paramSweep4.mat') % DBS data set
    case 4
        load('ecb43e_RHI_async_trial14.mat') % rubber hand illusion data set
    case 5
        load('3f2113_stim_12_52.mat') % stimulation spacing data set
end

%% 1/22/2018 - template subtraction
pre = 0.9; % started with 0.7
post = 1.8; % started with 0.7, then 1.1, then 1.5, then 1.8
[processedSig,templateDict_cell,template,start_inds,end_inds] = templateSubtract(dataInt,'type','average','fs',fs_data,'plotIt',0,'pre',pre,'post',post,'stimChans',stimChans);
%processedSig = templateSubtract(dataInt,fs_data,'plotIt',0);

numChans = size(dataInt,2);
[goods,goodVec] = goodChannel_extract('stimchans',stimChans,'numChans',numChans);

%
figure
hold on

chanInt = 2;
for i = 1:size(template{1},2)
    plot(template{chanInt}{i});
end
%
figure
hold on
for j = goodVec
    subplot(8,8,j)
    hold on
    for i = 1:size(templateDict_cell{j},2)
        timeVec = [0:size(templateDict_cell{j},1)-1];
        
        %plotBTLError(timeVec,templateDict_cell{j}(:,i),'CI');
        plot(timeVec,templateDict_cell{j}(:,i),'linewidth',2);
    end
    title(['Channel ' num2str(j)])
end

%% 1/22/2018 - linear interpolation

processedSig = interpolate_artifact(dataInt,'fs',fs_data,'plotIt',0,'type','pchip','stimchans',stimChans);
%%
avgResponse = mean(processedSig,3);
smallMultiples_responseTiming(avgResponse,t_epoch,'type1',stimChans,'type2',0,'average',1)
%%
processedSig_rms = rms_func(avgResponse(t_epoch<1000 & t_epoch>-100,:));
processedSig_var = var(avgResponse(t_epoch<1000 & t_epoch>-100,:));
figure
histogram(processedSig_rms)
[~,ind] = max(processedSig_rms);
ind
%% 1-29-2018 - template dictionary method

pre = 1; % started with 0.7
post = 2; % started with 0.7, then 1.1, then 1.5, then 1.8
[processedSig,templateDict_cell,template,start_inds,end_inds] = templateSubtract_dictionary(dataInt,'fs',fs_data,...
    'plotIt',0,'pre',pre,'post',post,'stimChans',stimChans);

chanIntList = [12 21 28 19 18 36 44 43 30 33 41 34];

for ind = chanIntList
    
    exampChan = mean(squeeze(processedSig(:,ind,:)),2);
    
    figure
    ax1 = subplot(2,1,1);
    plot(1e3*t_epoch,exampChan,'linewidth',2);
    xlim([-200 1000])
    ylim([-5e-5 5e-5])
    
    title(['Processed Signal - Channel ' num2str(ind)])
    clear exampChan
    
    
    ax2 = subplot(2,1,2);
    exampChan = mean(squeeze(dataInt(:,ind,:)),2);
    plot(1e3*t_epoch,exampChan,'linewidth',2);
    xlim([-200 1000])
    ylim([-5e-5 5e-5])
    xlabel('time (ms)')
    ylabel('Voltage (\muV)')
    title(['Raw Signal Average - Channel ' num2str(ind)])
    
    linkaxes([ax1,ax2],'xy')
    
    clear exampChan
end
%%=

figure
hold on
for j = goodVec
    subplot(8,8,j)
    hold on
    for i = 1:size(templateDict_cell{j},2)
        timeVec = [0:size(templateDict_cell{j},1)-1];
        
        %plotBTLError(timeVec,templateDict_cell{j}(:,i),'CI');
        plot(timeVec,templateDict_cell{j}(:,i),'linewidth',2);
    end
    title(['Channel ' num2str(j)])
end

%
sig = processedSig;
avgResponse = mean(sig,3);

stimChans = [20 29];

smallMultiples_responseTiming(avgResponse,t_epoch,'type1',stimChans,'type2',0,'average',1,'highlight_range',[0 800])

% plot the average dictionary templates
figure
for j = goodVec
    subplot(8,8,j)
    hold on
    for i = 1:size(template{j},2)
        timeVec = [0:size(template{j}{:,i},1)-1];
        
        plotBTLError(timeVec,template{j}{:,i},'CI',rand(1,3)');
        %plot(timeVec,templateDict_cell{j}(:,i),'linewidth',2);
    end
    title(['Channel ' num2str(j)])
end

%% - template trial method
stimChans = [1 9 20 24 29 32]; % 29 was bad too, 1 9 29 32 were the stim channels, 20 might also be bad
pre = 1; % started with 0.7, then 0.9, adjusted the window, move it to 0.4
post = 2; % started with 0.7, then 1.1, then 1.5, then 1.8
[processedSig_trial,templateDict_cell_trial,template_trial,start_inds,end_inds] = templateSubtract_trial(dataInt,'fs',fs_data,...
    'plotIt',0,'pre',pre,'post',post,'stimChans',stimChans);

numChans = size(dataInt,2);
[goods,goodVec] = goodChannel_extract('stimchans',stimChans,'numChans',numChans);

%
chanIntList = [12 21 28 19 18 36 44 43 30 33 41 34];

for ind = chanIntList
    
    exampChan = mean(squeeze(processedSig_trial(:,ind,:)),2);
    
    figure
    ax1 = subplot(2,1,1);
    plot(1e3*t_epoch,exampChan,'linewidth',2);
    xlim([-200 1000])
    ylim([-5e-5 5e-5])
    
    title(['Processed Signal - Channel ' num2str(ind)])
    clear exampChan
    
    
    ax2 = subplot(2,1,2);
    exampChan = mean(squeeze(dataInt(:,ind,:)),2);
    plot(1e3*t_epoch,exampChan,'linewidth',2);
    xlim([-200 1000])
    ylim([-5e-5 5e-5])
    xlabel('time (ms)')
    ylabel('Voltage (\muV)')
    title(['Raw Signal Average - Channel ' num2str(ind)])
    
    linkaxes([ax1,ax2],'xy')
    
    clear exampChan
end
%
sig = processedSig_trial;
avgResponse = mean(sig,3);

stimChans = [20 29];

smallMultiples_responseTiming(avgResponse,t_epoch,'type1',stimChans,'type2',0,'average',1)

% plot the average dictionary templates
figure
for j = goodVec
    subplot(8,8,j)
    hold on
    for i = 1:size(template_trial{j},2)
        timeVec = [0:size(template_trial{j}{:,i},1)-1];
        
        plotBTLError(timeVec,template_trial{j}{:,i},'CI',rand(1,3)');
        %plot(timeVec,templateDict_cell{j}(:,i),'linewidth',2);
    end
    title(['Channel ' num2str(j)])
end

%%
stimChans = [1 9 20 24 29 32]; % 29 was bad too, 1 9 29 32 were the stim channels, 20 might also be bad
pre = 0.3; % started with 0.7, then 0.9, adjusted the window, move it to 0.4
post = 0.8; % started with 0.7, then 1.1, then 1.5, then 1.8
[processedSig_v2,templateDict_cell_v2,template_v2] = templateSubtract_dictionary_iterative(processedSig,'fs',fs_data,'plotIt',0,...
    'pre',pre,'post',post,'stimChans',stimChans,'start_inds',start_inds,'end_inds',end_inds);


%%
chanIntList = [12 21 28 19 18 36 44 43 30 33 41 34];

for ind = chanIntList
    
    exampChan = mean(squeeze(processedSig_v2(:,ind,:)),2);
    
    figure
    ax1 = subplot(2,1,1);
    plot(1e3*t_epoch,exampChan,'linewidth',2);
    xlim([-200 1000])
    ylim([-5e-5 5e-5])
    
    title(['Processed Signal - Channel ' num2str(ind)])
    clear exampChan
    
    
    ax2 = subplot(2,1,2);
    exampChan = mean(squeeze(dataInt(:,ind,:)),2);
    plot(1e3*t_epoch,exampChan,'linewidth',2);
    xlim([-200 1000])
    ylim([-5e-5 5e-5])
    xlabel('time (ms)')
    ylabel('Voltage (\muV)')
    title(['Raw Signal Average - Channel ' num2str(ind)])
    
    linkaxes([ax1,ax2],'xy')
    
    clear exampChan
end

figure
hold on
for j = goodVec
    subplot(8,8,j)
    hold on
    for i = 1:size(templateDict_cell_v2{j},2)
        timeVec = [0:size(templateDict_cell_v2{j},1)-1];
        
        %plotBTLError(timeVec,templateDict_cell{j}(:,i),'CI');
        plot(timeVec,templateDict_cell_v2{j}(:,i),'linewidth',2);
    end
    title(['Channel ' num2str(j)])
end

%
sig = processedSig_v2;
avgResponse = mean(sig,3);

stimChans = [20 29];

smallMultiples_responseTiming(avgResponse,t_epoch,'type1',stimChans,'type2',0,'average',1)

% plot the average dictionary templates
figure
for j = goodVec
    subplot(8,8,j)
    hold on
    for i = 1:size(template_v2{j},2)
        timeVec = [0:size(template_v2{j}{:,i},1)-1];
        
        plotBTLError(timeVec,template_v2{j}{:,i},'CI',rand(1,3)');
        %plot(timeVec,templateDict_cell{j}(:,i),'linewidth',2);
    end
    title(['Channel ' num2str(j)])
end

%% PROCESS THE DATA
% process the wavelet using morlet process and PLV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trial by trial wavelet decomp, PLV

%%%%%% PLV
freq_range = [8 12];
[plv] = plvWrapper(processedSig,fs_data,freq_range,stimChans);

%%%%%%% wavelet
time_res = 0.050; % 50 ms bins

[powerout,f_morlet,t_morlet,~] = waveletWrapper(processedSig,fs_data,time_res,stimChans);
%
t_morlet = linspace(-pre_stim,post_stim,length(t_morlet))/1e3;
badChans = [1 9 20 24 29 32];
powerout(:,:,badChans,:) = 0;

% Visualize wavelets

% example wavelet decomp
%trialInt = 20; % which channel to check out
chanInt = 27;

%t_epoch = (-samps_pre_stim:samps_post_stim-1)/fs_data;

response = buttonLocs{condInt}/(2*fs_data);
cmap=flipud(cbrewer('div', 'RdBu', 13));
colormap(cmap)

powerout_norm = normalize_spectrogram(powerout,t_morlet);


for i = 1:size(powerout,4)
    
    totalFig = figure;
    totalFig.Units = 'inches';
    totalFig.Position = [12.1806 3.4931 6.0833 7.8056];
    subplot(3,1,1);
    imagesc(1e3*t_morlet,f_morlet,powerout_norm(:,:,chanInt,i));
    cmap=flipud(cbrewer('div', 'RdBu', 13));
    colormap(cmap)
    axis xy;
    xlabel('time (ms)');
    ylabel('frequency (Hz)');
    title(['Wavelet decomposition Channel ' num2str(chanInt) ' Trial ' num2str(i)]);
    vline(stimTime(i),'m','stim')
    vline(1e3*response(i),'g','response')
    %xlim([-200 1000]);
    %colorbar()
    caxis([-3 3])
    set(gca,'fontsize',14)
    
    
    
    %figure;
    h1 = subplot(3,1,2);
    plot(1e3*t_epoch,1e6*processedSig(:,chanInt,i))
    vline(stimTime(i),'m','stim')
    xlabel('time (ms)');
    ylabel('microvolts')
    title(['Processed Channel ' num2str(chanInt) ' Trial ' num2str(i)]);
    vline(1e3*response(i),'g','response')
    ylims = [-(max(abs(1e6*processedSig(:,chanInt,i))) + 10) (max(abs(1e6*processedSig(:,chanInt,i))) + 10)];
    ylim(ylims);
    ylim_h1 = ylims;
    %xlim([-200 1000]);
    set(gca,'fontsize',14)
    
    
    h2 = subplot(3,1,3);
    plot(1e3*t_epoch,1e6*dataInt(:,chanInt,i))
    vline(stimTime(i),'m','stim')
    xlabel('time (ms)');
    ylabel('microvolts')
    title(['Raw Channel ' num2str(chanInt) ' Trial ' num2str(i)]);
    vline(1e3*response(i),'g','response')
    ylim(ylim_h1);
    %xlim([-200 1000]);
    set(gca,'fontsize',14);
    
    
    
    linkaxes([h1,h2],'xy');
    
end
% normalized spectogram

powerout_norm = normalize_spectrogram(powerout,t_morlet);

avg_power_norm = mean(powerout_norm,4);
stimChans = badChans;
smallMultiples_responseTiming_spectrogram(avg_power_norm,t_morlet,f_morlet,'type1',stimChans,'type2',0,'average',1);


