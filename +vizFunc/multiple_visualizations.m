function multiple_visualizations(processedSig,rawSig,varargin)
% Usage:  [processedSig,startInds,endInds] = interpolate_artifact(rawSig,varargin)
%
% Function to visualize
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get inputs
p = inputParser;

validData = @(x) isnumeric(x);
addRequired(p,'processedSig',validData);
addRequired(p,'rawSig',validData);

addParameter(p,'type','linear',@isstr);

addParameter(p,'xlims',[-100 1000],@isnumeric);
addParameter(p,'ylims',[-6 6],@isnumeric);
addParameter(p,'trainDuration',[0 400],@isnumeric);
addParameter(p,'tEpoch',0.2,@isnumeric);
addParameter(p,'stimChans',[],@isnumeric);
addParameter(p,'bads',[],@isnumeric);
addParameter(p,'chanIntList',[1,2,3],@isnumeric);
addParameter(p,'fs',12207,@isnumeric);
addParameter(p,'templateTrial',{},@iscell)
addParameter(p,'templateDictCell',{},@iscell);
addParameter(p,'modePlot','avg',@isstr);

p.parse(processedSig,rawSig,varargin{:});

rawSig = p.Results.rawSig;

type = p.Results.type;

stimChans = p.Results.stimChans;
bads = p.Results.bads;
fs = p.Results.fs;

templateTrial = p.Results.templateTrial;
templateDictCell = p.Results.templateDictCell;
trainDuration = p.Results.trainDuration;
xlims = p.Results.xlims;
ylims = p.Results.ylims;
tEpoch = p.Results.tEpoch;
chanIntList = p.Results.chanIntList;
modePlot = p.Results.modePlot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
numChans = size(rawSig,2);
[goods,goodVec] = helpFunc.good_channel_extract('bads',bads,'stimchans',stimChans,'numChans',numChans);
p = vizFunc.numSubplots(numChans);

if (strcmp(type,'dictionary') || strcmp(type,'trial') || strcmp(type,'average')) && (~isempty(templateDictCell) || ~isempty(templateTrial))
    
    if exist('templateTrial','var')
        figure
        for j = goodVec
            subplot(p(1),p(2),j)
            hold on
            for i = 1:size(templateTrial{j},2)
                timeVec = 1e3*[0:size(templateTrial{j}{:,i},1)-1]/fs;
                
                vizFunc.plot_error(timeVec,1e3*templateTrial{j}{:,i},'CI',rand(1,3));
            end
            title(['Channel ' num2str(j)])
        end
        xlabel('Time (ms')
        ylabel('Voltage (mV)')
    end
    
    % plot the average template dictionary if using a dictionary method
    if strcmp(type,'dictionary') && exist('templateDictCell','var')
        figure
        
        hold on
        for j = goodVec
            subplot(p(1),p(2),j)
            hold on
            for i = 1:size(templateDictCell{j},2)
                timeVec = 1e3*[0:size(templateDictCell{j},1)-1]/fs;
                
                plot(timeVec,1e3*templateDictCell{j}(:,i),'linewidth',2);
            end
            title(['Channel ' num2str(j)])
        end
        xlabel('Time (ms')
        ylabel('Voltage (mV)')
    end
end

avgResponse = mean(processedSig,3);
avgRaw = mean(rawSig,3);

vizFunc.small_multiples_time_series(processedSig(:,:,:),tEpoch,'type1',stimChans,'type2',0,'xlims',xlims,'modePlot',modePlot,'highlightRange',trainDuration)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make a list of the channels of interest to visualize one at a time
%%
for ind = chanIntList
    
    exampChan = mean(squeeze(processedSig(:,ind,:)),2);
    
    figure
    ax1 = subplot(2,1,1);
    plot(1e3*tEpoch,1e3*exampChan,'linewidth',2,'color',[204 85 0]/255);
    xlim(xlims)
    ylim(ylims)
    set(gca,'fontsize',18)
    title(['Recovered Signal - Channel ' num2str(ind)])
    clear exampChan
    
    
    ax2 = subplot(2,1,2);
    exampChan = mean(squeeze(rawSig(:,ind,:)),2);
    plot(1e3*tEpoch,1e3*exampChan,'linewidth',2,'color','k');
    xlim(xlims)
    ylim(ylims)
    xlabel('Time (ms)')
    ylabel('Voltage (mV)')
    title(['Raw Signal - Channel ' num2str(ind)])
    set(gca,'fontsize',18)
    
    linkaxes([ax1,ax2],'xy')
    
    clear exampChan
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% look at the FFT difference
[f,P1] = helpFunc.fourier_transform_calc(fs,avgResponse);
[fRaw,P1Raw] = helpFunc.fourier_transform_calc(fs,avgRaw);

vizFunc.small_multiples_fourier(P1Raw,fRaw,'type1',stimChans,'type2',0)
legend('raw')
vizFunc.small_multiples_fourier(P1,f,'type1',stimChans,'type2',0,'newfig',0)
legend('processed')

vizFunc.small_multiples_fourier(P1Raw,fRaw,'type1',stimChans,'type2',0,'plotLog',1)
legend('raw')
vizFunc.small_multiples_fourier(P1,f,'type1',stimChans,'type2',0,'newfig',0,'plotLog',1)
legend('processed')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% look at the RMS reduction within particular section
processedSigRms = helpFunc.rms_func(avgResponse(tEpoch<1000 & tEpoch>-100,:));
rawSigRms = helpFunc.rms_func(avgRaw(tEpoch<1000 & tEpoch>-100,:));

% look at decibel reduction for each channel

rmsDb = 20*log10(processedSigRms./rawSigRms);

figure
numBins = 10;
histogram(rmsDb,numBins)
title('RMS Reduction in Decibels')
xlabel('magnitude of decibel decrease')
ylabel('count')
set(gca,'fontsize',14)

end