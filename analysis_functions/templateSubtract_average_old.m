function [processedSig,template] = templateSubtract(raw_sig,varargin)
%USAGE:
% This function will perform a template subtraction scheme for artifacts on
% a trial by trial, channel by channel basis
%
% raw_sig = samples x channels x trials
% pre = the number of ms before which the stimulation pulse onset as
% detected by a thresholding method should still be considered as artifact
% post = the number of ms after which the stimulation pulse onset as
% detected by a thresholding method should still be considered as artifact
% fs = sampling rate (Hz)
% plotIt = plot intermediate steps if true

% default parameters
fs = 1.2207e04;
plotIt = 0;
type = 'linear';
pre = 0.4096;
post = 0.4096;

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
        case 'stimchans'
            stimChans = varargin{i+1};
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

% take diff of signal to find onset of stimulation train
diff_sig = permute(cat(3,zeros(size(raw_sig,2), size(raw_sig,3)),permute(diff(raw_sig),[2 3 1])),[3 1 2]);

% find channel that has the max signal, and use this for subsequent
% analysis
[~,chanMax] = (max(max(diff_sig)));
chanMax = chanMax(1);
lengthMax = 0; % length vector to build up the dictionary of templates later
template = {};

for trial = 1:size(raw_sig,3)
    
    % find onset
    
%     inds = find(abs(zscore(diff_sig(:,chanMax,trial)))>0.5);
%     %inds = find(abs(raw_sig(:,chanMax,trial))>2e-4);
%     diff_bt_inds = [diff([0; inds])'];
%     [~,inds_onset] = find(abs(zscore(diff_bt_inds))>0.5);
%     if length(inds_onset == 2)
%         start_inds = inds(inds_onset(1))-2;
%         end_inds = inds(inds_onset(2)-1)+2;
%     else
%         start_inds = [inds(1)-presamps; inds(inds_onset+1)-presamps];
%         end_inds = [ inds(inds_onset)+postsamps; inds(end)+postsamps];
%     end
    
    
    inds = find(abs(zscore(diff_sig(:,chanMax,trial)))>2);
    %inds = find(diff_sig(:,chanMax,trial)>2e-4);
    diff_bt_inds = [diff(inds)'];
    [~,inds_onset] = find(abs(zscore(diff_bt_inds))>2);
    %[~,inds_onset] = find(diff_bt_inds>5);
    start_inds{trial} = [inds(1)-presamps; inds(inds_onset+1)-presamps];
    end_inds{trial} = start_inds{trial}+postsamps+17;
    %end_inds = [ inds(inds_onset)+postsamps; inds(end)+postsamps];
    
    
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
        % get single trial
        raw_sig_temp = raw_sig(:,chan,trial);
        % default time to average
        default_time_average = end_inds{trial}(1)-start_inds{trial}(1)+1;
                lengthMax = max(default_time_average,lengthMax);

        avg_signal = zeros((end_inds(1)-start_inds(1)+1),length(start_inds));
        for sts = 1:length(start_inds)
            win = start_inds{trial}(sts):start_inds{trial}(sts)+default_time_average-1;
            % avg_signal(:,sts) = raw_sig_temp(win);
            avg_signal(:,sts) = raw_sig_temp(win) - mean(raw_sig_temp(start_inds{trial}(sts):start_inds{trial}(sts)+presamps-8));% take off pre period
        end
        
        template{chan}{trial} = avg_signal;
    
    end
end

for chan = goodVec
    templateArray = [];
    
    for trial = 1:size(raw_sig,3)
        artifacts = template{chan}{trial};
        
        if size(artifacts,1) < lengthMax
            amntPad = lengthMax - size(artifacts,1);
            artifacts = padarray(artifacts,[amntPad 0],'post');
        end
        
        templateArray = [templateArray artifacts];
        
    end
    
    avg_signal_mean = mean(templateArray,2);
    
    for trial = 1:size(raw_sig,3)
        raw_sig_temp = raw_sig(:,chan,trial);

        for sts = 1:length(start_inds)
            win = start_inds{trial}(sts):start_inds{trial}(sts)+default_time_average-1;
            raw_sig_temp(win) = raw_sig_temp(win) - avg_signal_mean;
        end
        
        processedSig(:,chan,trial) = raw_sig_temp;
    end
    
end


end