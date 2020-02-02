function [startInds,endInds] = get_artifact_indices(rawSig,varargin)

% Usage:  [startInds,endInds] = get_artifact_indices(rawSig,varargin)
%
% This function will extract the indices to begin and end each artifact
% selection period on a channel and trial basis. The channel with the
% largest artifact is used to select the approximate beginning of the
% artifacts across all other channels.
%
% Arguments:
%   Required:
%   rawSig - samples x channels x trials
%
%   Optional:
%useFixedEnd - use a fixed end distance (1), or dynamically calculate the
%              offset of each stimulus pulse
%      fixedDistance - the maximum distance in ms to either look beyond
%        pre - the number of ms before which the stimulation pulse onset as
%              detected by a thresholding method should still be considered
%              as artifact
%       post - the number of ms after which the stimulation pulse onset as
%              detected by a thresholding method should still be
%              considered as artifact
%  preInterp - the number of ms before the stimulation which to consider an
%              interpolation scheme on. Does not apply to the linear case
% postInterp - the number of ms before the stimulation which to consider an
%              interpolation scheme on. Does not apply to the linear case
%          fs - sampling rate (Hz)
%      plotIt - plot intermediate steps if true
%     goodVec - vector (e.g. [1 2 4 10]) of the channel indices to use. If
%              there are bad channels or stimulation channels for instance,
%              they should be excluded from goodVec
%onsetThreshold - This value is used as absolute valued z-score threshold
%               to determine the onset of artifacts within a train. The differentiated
%               smoothed signal is used to determine artifact onset. This is also used in
%               determining how many stimulation pulses are within a train, by ensuring
%               that beginning of each artifact is within a separate artifact pulse.
%threshVoltageCut - This is used to help determine the end of each
%               individual artifact pulse dynamically. More specifically, this is a
%               percentile value, against which the absolute valued, z-scored smoothed raw
%               signal is compared to find the last value which exceeds the specified
%               percentile voltage value. Higher values of this (i.e. closer to 100)
%               result in a shorter duration artifact, while lower values result in a
%               longer calculated duration of artifact. This parameter therefore should
%               likely be set higher for more transient artifacts and lower for longer
%               artifacts.
%threshDiffCut - This is used to help determine the end of each individual
%               artifact pulse dynamically. More specifically, this is a percentile value,
%               against which the absolute valued, z-scored differentiated smoothed raw
%               signal is compared to find the last value which exceeds the specified
%               percentile voltage value. Higher values of this (i.e. closer to 100)
%               result in a shorter duration artifact, while lower values result in a
%               longer calculated duration of artifact. This parameter therefore should
%               likely be set higher for more transient artifacts and lower for longer
%               artifacts.
%
%
% Returns:
%      startInds - cell array of the start indices each artifact for each
%      channel and trial - startInds{trial}{channel}
%       endsInds - cell array of the end indices of each artifact for each
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
p = inputParser;

validData = @(x) isnumeric(x);
addRequired(p,'rawSig',validData);

addParameter(p,'useFixedEnd',0,@(x) x==0 || x ==1);
addParameter(p,'fixedDistance',2,@isnumeric);

addParameter(p,'pre',0.4096,@isnumeric);
addParameter(p,'plotIt',0,@(x) x==0 || x ==1);
addParameter(p,'post',0.4096,@isnumeric);
addParameter(p,'fs',12207,@isnumeric);
addParameter(p,'goodVec',[1:64],@isnumeric);
addParameter(p,'chanInt',1,@isnumeric);
addParameter(p,'minDuration',0,@isnumeric);

addParameter(p,'onsetThreshold',1.5,@isnumeric);

addParameter(p,'threshVoltageCut',75,@isnumeric);
addParameter(p,'threshDiffCut',75,@isnumeric);

p.parse(rawSig,varargin{:});

rawSig = p.Results.rawSig;
plotIt = p.Results.plotIt;
useFixedEnd = p.Results.useFixedEnd;
pre = p.Results.pre;
post = p.Results.post;
fixedDistance = p.Results.fixedDistance;
fs = p.Results.fs;
goodVec = p.Results.goodVec;
chanInt = p.Results.chanInt;
minDuration = p.Results.minDuration;
threshVoltageCut = p.Results.threshVoltageCut;
threshDiffCut = p.Results.threshDiffCut;
onsetThreshold = p.Results.onsetThreshold;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

presamps = round(pre/1e3 * fs); % pre time in sec
postsamps = round(post/1e3 * fs); %
minDuration = round(minDuration/1e3 * fs);
fixedDistanceSamps = round(fixedDistance/1e3 * fs);
defaultWinAverage = fixedDistanceSamps ;

% take diff of signal to find onset of stimulation train
order = 3;
framelen = 7;
fprintf(['-------Smoothing data with Savitsky-Golay Filter-------- \n'])

rawSigFilt = rawSig;

for ind2 = 1:size(rawSigFilt,2)
    for ind3 = 1:size(rawSigFilt,3)
        rawSigFilt(:,ind2,ind3) = savitskyGolay.sgolayfilt_complete(squeeze(rawSig(:,ind2,ind3)),order,framelen);
        
    end
end

diffSig = permute(cat(3,zeros(size(rawSig,2), size(rawSig,3)),permute(diff(rawSigFilt),[2 3 1])),[3 1 2]);

fprintf(['-------Done smoothing and differentiating-------- \n'])


% find channel that has the max signal, and use this for subsequent
% analysis
[~,chanMax] = (max(max(diffSig(:,goodVec,:))));
chanMax = chanMax(1);
chanMax = goodVec(chanMax);
lengthMax_vec = []; % length vector to build up the dictionary of templates later
fprintf(['-------Getting artifact indices-------- \n'])

if ~exist('chanInt','var')
    chanInt = chanMax;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pctl = @(v,p) interp1(linspace(0.5/length(v), 1-0.5/length(v), length(v))', sort(v), p*0.01, 'spline');

timeSampsExtend = 2*fs/1000;% time_ms

for trial = 1:size(rawSig,3)
    
    inds = find(abs(zscore(diffSig(:,chanMax,trial)))>onsetThreshold); 
    diffBtInds = [diff(inds)'];
    [~,indsOnset] = find(abs(zscore(diffBtInds))>onsetThreshold);
    
    for chan = goodVec
        
        startInds{trial}{chan} = [inds(1)-presamps; inds(indsOnset+1)-presamps]';
        
        if useFixedEnd
            endInds{trial}{chan} = startInds{trial}{chan}+fixedDistanceSamps;
        else
            
            for idx = 1:length(startInds{trial}{chan})
                
                win = startInds{trial}{chan}(idx):startInds{trial}{chan}(idx)+defaultWinAverage; % get window that you know has the end of the stim pulse
                signal = rawSigFilt(win,chan,trial);
                diffSignal = diffSig(win,chan,trial);
          
                absZSig = abs(zscore(signal));
                absZDiffSig = abs(zscore(diffSignal));

                threshSig = pctl(absZSig,threshVoltageCut); 
                threshDiff = pctl(absZDiffSig,threshDiffCut); 
                
                % look past minimum start time
                last = presamps+minDuration+find(absZSig(presamps+minDuration:end)>threshSig,1,'last'); 
                last2 = presamps+minDuration+find(absZDiffSig(presamps+minDuration:end)>threshDiff,1,'last')+1; 
                ct = max(last, last2);
                
                if isempty(ct)
                    ct = last;
                    if isempty(last)
                        ct = last2;
                        if isempty(last2)
                            ct = postsamps;
                        end
                    end
                end
                
                endInds{trial}{chan}(idx) = ct + startInds{trial}{chan}(idx) + postsamps;
                
            end
        end
        
        if plotIt
            %%
            tDiff = 1e3*(0:length(diffSig(:,chanInt,trial))-1)/fs;
            tRaw = 1e3*(0:length(rawSig(:,chanInt,trial))-1)/fs;
            figure
            p1 =   subplot(2,1,1);
            plot(tDiff,abs(zscore(diffSig(:,chanInt,trial))),'linewidth',2,'color','k');
            h1 = vline(1e3*startInds{trial}{chanInt}/fs);
            h2 = vline(1e3*endInds{trial}{chanInt}/fs,'b:');
            ylabel('Z-score')
            xlabel('Time (ms)')
            title('Absolute Value Z-Scored Differentiated Smoothed Signal')
            set(gca,'fontsize',18)
            legend([h1(1) h2(1)],{'beginning','end'});
            
            p2 =   subplot(2,1,2);
            plot(tRaw,abs(zscore(rawSig(:,chanInt,trial))),'linewidth',2,'color','k');
            vline(1e3*startInds{trial}{chanInt}/fs)
            vline(1e3*endInds{trial}{chanInt}/fs,'b:')
            ylabel('Z-score')
            xlabel('Time (ms)')
            title('Absolute Value Z-Scored Raw Signal')
            set(gca,'fontsize',18)
            
            linkaxes([p1 p2],'xy')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            figure
            p3=    subplot(2,1,2);
            plot(tDiff,diffSig(:,chanInt,trial),'linewidth',2,'color','k');
            h1 = vline(1e3*startInds{trial}{chanInt}/fs);
            h2 = vline(1e3*endInds{trial}{chanInt}/fs,'b:');
            ylabel('Voltage (mV)')
            xlabel('Time (ms)')
            title('Differentiated Smoothed Signal')
            set(gca,'fontsize',18)
            legend([h1(1) h2(1)],{'beginning','end'});
            
            p4 =   subplot(2,1,1);
            plot(tRaw,rawSig(:,chanInt,trial),'linewidth',2,'color','k');
            vline(1e3*startInds{trial}{chanInt}/fs)
            vline(1e3*endInds{trial}{chanInt}/fs,'b:')
            ylabel('Voltage (mV)')
            xlabel('Time (ms)')
            title('Raw Signal')
            set(gca,'fontsize',18)
            linkaxes([p3 p4],'xy')
            
        end
        
    end
    fprintf(['-------Finished getting artifacts - Trial ' num2str(trial) '-------- \n'])
    
end

end