function [processedSig] = interpolate_artifact(raw_sig,varargin)
%USAGE:
% This function will perform an interpolation scheme for artifacts on a
% trial by trial, channel by channel basis, implementing either a linear
% interpolation scheme, or a pchip interpolation scheme
%
% raw_sig = samples x channels x trials
% pre = the number of ms before which the stimulation pulse onset as
% detected by a thresholding method should still be considered as artifact
% post = the number of ms after which the stimulation pulse onset as
% detected by a thresholding method should still be considered as artifact
% pre_interp = the number of ms before the stimulation which to consider an
% interpolation scheme on. Does not apply to the linear case
% post_interp = the number of ms before the stimulation which to consider an
% interpolation scheme on. Does not apply to the linear case
% fs = sampling rate (Hz)
% plotIt = plot intermediate steps if true

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check some basic data requirements
if nargin == 0
  error ('You must supply data');
end

if length (size (raw_sig)) > 3
  error ('Input data can not have more than three dimensions.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default parameters
fs = 1.2207e04;
plotIt = 0;
type = 'linear';
pre = 0.4096;
post = 0.4096;
pre_interp = 0.2; %0.5734;
post_interp = 0.2;
%pre_interp = 0.9;
%post_interp = 0.9;
stimChans = [];
useFixedEnd = 1;
fixed_distance = 1;

% define matrix of zeros
processedSig = zeros(size(raw_sig));

for i=1:2:(length(varargin)-1)
    
    switch lower(varargin{i})
        case 'plotit'
            plotIt = varargin{i+1};
        case 'fs'
            fs = varargin{i+1};
        case 'type'
            type = varargin{i+1};
        case 'pre'
            pre = varargin{i+1};
        case 'post'
            post = varargin{i+1};
        case 'pre_interp'
            pre_interp = varargin{i+1};
        case 'post_interp'
            post_interp = varargin{i+1};
        case 'stimchans'
            stimChans = varargin{i+1};
        case 'usefixedend'
            useFixedEnd = varargin{i+1};
        case 'fixed_distance'
            fixed_distance = varargin{i+1};
    end
end

bads = [];
badTotal = [stimChans; bads];

% total channels
numChans = size(raw_sig,2);
% make logical good channels matrix to index
goods = zeros(numChans,1);

channelsOfInt = 1:numChans;

goods(channelsOfInt) = 1;
% set the goods matrix to be zero where bad channels are
goods(badTotal) = 0;
% make it logical
goods = logical(goods);

%raw_sig = raw_sig(:,goods,:);

[~,goodVec] = find(goods'>0);

%% preprocess eco
presamps = round(pre/1e3 * fs); % pre time in sec

postsamps = round(post/1e3 * fs); %

pre_interp_samps = round(pre_interp/1e3 * fs);

post_interp_samps = round(post_interp/1e3 * fs);

fixed_distance_samps = round(fixed_distance/1e3 * fs); 
default_win_average = fixed_distance_samps ; %end_inds{trial}(1)-start_inds{trial}(1)+1;

% take diff of signal to find onset of stimulation train
diff_sig = permute(cat(3,zeros(size(raw_sig,2), size(raw_sig,3)),permute(diff(raw_sig),[2 3 1])),[3 1 2]);

% find channel that has the max signal, and use this for subsequent
% analysis
[~,chanMax] = (max(max(diff_sig(:,goodVec,:))));
chanMax = chanMax(1);

fprintf(['-------Interpolation-------- \n'])
fprintf(['-------' type '-------- \n'])

for trial = 1:size(raw_sig,3)

    inds = find(abs(zscore(diff_sig(:,chanMax,trial)))>2);
    diff_bt_inds = [diff(inds)'];
    [~,inds_onset] = find(abs(zscore(diff_bt_inds))>2);
    start_inds{trial} = [inds(1)-presamps; inds(inds_onset+1)-presamps];
    
    if useFixedEnd
        end_inds{trial} = start_inds{trial}+fixed_distance_samps; % 17 to start
    else
        
        for idx = 1:length(start_inds{trial})
            
            win = start_inds{trial}(idx):start_inds{trial}(idx)+default_win_average; % get window that you know has the end of the stim pulse 
            signal = raw_sig(win,chanMax,trial);
            diff_signal = diff_sig(win,chanMax,trial);

            last = find(abs(zscore(signal))>0.2,1,'last');
            last2 = find(abs(zscore(diff_signal))>5e-3,1,'last')+1;
            
            if (isempty(last2))
                if (isempty(last))
                    error ('something is wrong with signal');
                else
                    ct = last;
                end
            else
                if (isempty(last))
                    ct = last2;
                else
                    ct = max(last, last2);
                end
            end
             end_inds{trial}(idx) = ct + start_inds{trial}(idx) + postsamps;

        end
    end
  
    if plotIt
        figure
        plot(diff_sig(:,chanMax,trial))
        vline(start_inds{trial})
        vline(end_inds{trial},'g')
        figure
        plot(raw_sig(:,chanMax,trial))
        vline(start_inds{trial})
        vline(end_inds{trial},'g')
    end
    
    for chan = goodVec
        raw_sig_temp = raw_sig(:,chan,trial);
        for sts = 1:length(start_inds{trial})
            win = start_inds{trial}(sts):end_inds{trial}(sts);
            switch type
                case 'linear'
                    raw_sig_temp(win) = interp1([start_inds{trial}(sts)-1 end_inds{trial}(sts)+1],...
                        raw_sig_temp([start_inds{trial}(sts)-1 end_inds{trial}(sts)+1]), start_inds{trial}(sts):end_inds{trial}(sts));
                case 'pchip'
                    raw_sig_temp(win) = interp1([start_inds{trial}(sts)-pre_interp_samps:start_inds{trial}(sts)-1 end_inds{trial}(sts):end_inds{trial}(sts)+post_interp_samps],...
                        raw_sig_temp([start_inds{trial}(sts)-pre_interp_samps:start_inds{trial}(sts)-1 end_inds{trial}(sts):end_inds{trial}(sts)+post_interp_samps]),...
                        start_inds{trial}(sts):end_inds{trial}(sts),'pchip');
            end
            
        end
        if plotIt && (trial == 10 || trial == 1000)
            figure
            plot(raw_sig_temp,'linewidth',2)
            hold on
            plot(raw_sig(:,chan,trial),'linewidth',2)
            
            vline(start_inds{trial})
            vline(end_inds{trial},'g')
        end
        processedSig(:,chan,trial) = raw_sig_temp;
        
    end
    
end

fprintf(['-------Finished-------- \n \n'])

end