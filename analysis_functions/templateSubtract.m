function [processedSig,templateArray_cell_output,template,start_inds,end_inds] = templateSubtract(raw_sig,varargin)
%USAGE:
% This function will perform a template subtraction scheme for artifacts on
% a trial by trial, channel by channel basis. This function will build up
% a dictionary of artifacts, best match the template, and
%
% raw_sig = samples x channels x trials
% pre = the number of ms before which the stimulation pulse onset as
% detected by a thresholding method should still be considered as artifact
% post = the number of ms after which the stimulation pulse onset as
% detected by a thresholding method should still be considered as artifact
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
% default parameters in case user does not supply any
fs = 1.2207e04;
plotIt = 1;
pre = 0.4096;
post = 0.4096;
stimChans = [];
bads = [];
distance_metric_dbscan = 'eucl';
distance_metric_sigMatch = 'eucl';
useFixedEnd = 1;

% define matrix of zeros
processedSig = zeros(size(raw_sig));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use input variables
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
        case 'bads'
            bads = varargin{i+1};
            
            % options are 'corr', 'eucl', 'cosine'
        case 'distance_metric_dbscan'
            distance_metric_dbscan = varargin{i+1};
        case 'distance_metric_sigmatch'
            distance_metric_sigMatch = varargin{i+1};
        case 'usefixedend'
            useFixedEnd = varargin{i+1};
        case 'fixed_distance'
            fixed_distance = varargin{i+1};
            
    end
end

% make a vector of the good channels to process
numChans = size(raw_sig,2);
[goods,goodVec] = goodChannel_extract('bads',bads,'stimchans',stimChans,'numChans',numChans);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preprocess eco
presamps = round(pre/1e3 * fs); % pre time in sec

postsamps = round(post/1e3 * fs); %

fixed_distance_samps = round(fixed_distance/1e3 * fs);

default_win_average = fixed_distance_samps ; %end_inds{trial}(1)-start_inds{trial}(1)+1;

% take diff of signal to find onset of stimulation train
diff_sig = permute(cat(3,zeros(size(raw_sig,2), size(raw_sig,3)),permute(diff(raw_sig),[2 3 1])),[3 1 2]);

% find channel that has the max signal, and use this for subsequent
% analysis
[~,chanMax] = (max(max(diff_sig(:,goodVec,:))));
chanMax = chanMax(1);
lengthMax_vec = []; % length vector to build up the dictionary of templates later
fprintf(['-------Templates-------- \n'])

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
        % get single trial
        raw_sig_temp = raw_sig(:,chan,trial);
        
        avg_signal = {};
        lengthMax_vec_temp = [];
        for sts = 1:length(start_inds{trial})
            win = start_inds{trial}(sts):end_inds{trial}(sts);
            lengthMax_temp = length(win);
             avg_signal{sts} = raw_sig_temp(win) - mean(raw_sig_temp(start_inds{trial}(sts):start_inds{trial}(sts)+presamps-8));% take off average of first 3 samples
            
            lengthMax_vec_temp = [lengthMax_vec_temp lengthMax_temp];
            
        end
        
        lengthMax_vec = [lengthMax_vec max(lengthMax_vec_temp)];
        template_cell{chan}{trial} = avg_signal;
        
        if plotIt && (trial == 10 || trial == 1000)
            figure
            plot(raw_sig_temp,'linewidth',2)
            vline(start_inds{trial})
            vline(end_inds{trial},'g')
        end
        
    end
    
end
lengthMax = max(lengthMax_vec);
fprintf(['-------Finished getting artifacts-------- \n'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get templates all same length
for chan = goodVec
    templateArray = [];
    templateArray_extracted = [];
    
    for trial = 1:size(raw_sig,3)
        artifacts_cell = template_cell{chan}{trial};
        artifacts_mat = [];
        for sts = 1:length(start_inds{trial})
            
            artifacts_trial = artifacts_cell{sts};
            
            if size(artifacts_trial,1) < lengthMax
                amntPad = lengthMax - size(artifacts_trial,1);
                artifacts_pad = padarray(artifacts_trial,amntPad,nan,'post');
            else
                artifacts_pad = artifacts_trial;
            end
            
            artifacts_mat(:,sts) = artifacts_pad;
            
        end
        template{chan}{trial} = artifacts_mat;
        templateArray = [templateArray artifacts_mat];
        templateArray_cell{chan} = templateArray;
    end
end

fprintf(['-------Finished making artifacts the same length-------- \n'])

%% build up dictionary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch type
    case 'dictionary'
        templateArray_cell_output = {};
        
        fprintf(['-------Dictionary-------- \n'])
        
        for chan = GoodVec
            
            templateArray = templateArray_cell{chan};         
            if max(templateArray(:)) < 3e-4
                distanceDBscan = 0.9;
            else
                distanceDBscan = 0.95;
            end
            
            [c,ptsC,centres] = dbscan(templateArray,distanceDBscan,1,distance_metric_dbscan);
            
            vectorUniq = unique(ptsC);
            if plotIt
                figure
                hold on
                vectorUniq = unique(ptsC)';
                for i = 1:length(vectorUniq)
                    plot(mean(templateArray(:,ptsC==i),2));
                end
            end
            
            
            for i = 1:length(vectorUniq)
                meanTempArray = mean(templateArray(:,ptsC==i),2);
                templateArray_extracted = [templateArray_extracted (meanTempArray )]; %no subtraction
            end
            
            templateArray_cell_output{chan} = templateArray_extracted;
            
            fprintf(['-------Artifact Channel ' num2str(chan) ' -------- \n'])
            
            
        end
        
        fprintf(['-------Finished clustering artifacts-------- \n'])
        
        %
        % now do the template subtraction
        for trial = 1:size(raw_sig,3)
            
            for chan = goodVec
                
                raw_sig_temp = raw_sig(:,chan,trial);
                
                for sts = 1:length(start_inds{trial})
                    win = start_inds{trial}(sts):end_inds{trial}(sts);
                    extracted_sig = raw_sig_temp(win);
                    extracted_sig = extracted_sig - extracted_sig(1);
                    
                    % find best artifact
                    
                    templates = templateArray_cell{chan};
                    
                    % add on the trial one
                    templates = [templates mean(template{chan}{trial},2)];
                    
                    % get them to be the same length
                    templates = templates(1:length(extracted_sig),:);
                    
                    switch distance_metric_sigMatch
                        case 'correlation'
                            % correlation
                            correlations = corr(extracted_sig,templates);
                            [~,index] = max((correlations));
                            
                        case 'cosineSim'
                            % - cosine similarity
                            denominator = sqrt(sum(extracted_sig.*extracted_sig)).*sqrt(sum(templates.*templates));
                            numerator = (extracted_sig'*templates);
                            correlations = numerator./denominator;
                            [~,index] = max(abs(correlations));
                            
                        case 'eucl'
                            % distance
                            v = templates - repmat(extracted_sig,1,size(templates,2));
                            distance = sum(v.*v);
                            [~,index] = min(distance);
                    end
                    
                    template_subtract = templates(:,index);
                    
                    raw_sig_temp(win) = raw_sig_temp(win) - template_subtract;
                end
                
                if plotIt && (trial == 10 || trial == 1000)
                    figure
                    plot(raw_sig_temp,'linewidth',2)
                    hold on
                    plot(raw_sig(:,chan,trial),'linewidth',2)
                    
                    vline(start_inds{trial})
                    vline(end_inds{trial},'g')
                    
                    figure
                    plot(extracted_sig)
                    hold on
                    plot(template_subtract)
                    plot(extracted_sig-template_subtract)
                    legend('extracted','template','subtracted');
                end
                
                processedSig(:,chan,trial) = raw_sig_temp;
                
            end
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'average'
        templateArray_cell_output = {};
        
        fprintf(['-------Average Template-------- \n'])
        
        for chan = goodVec
            templateArray = templateArray_cell{chan};
            
            avg_signal_mean = nanmean(templateArray,2);
            templateArray_cell{chan} = avg_signal_mean;
            
            for trial = 1:size(raw_sig,3)
                raw_sig_temp = raw_sig(:,chan,trial);
                
                for sts = 1:length(start_inds{trial})
                    win = start_inds{trial}(sts):end_inds{trial}(sts);
                    raw_sig_temp(win) = raw_sig_temp(win) - avg_signal_mean(1:length(win));
                end
                
                processedSig(:,chan,trial) = raw_sig_temp;
            end
            
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'trial'
        fprintf(['-------Trial Based Template-------- \n'])
        templateArray_cell_output = {};
        
        for trial = 1:size(raw_sig,3)
            for chan = goodVec
                raw_sig_temp = raw_sig(:,chan,trial);
                artifacts = template{chan}{trial};
                
                avg_trial = nanmean(artifacts,2);
                for sts = 1:length(start_inds{trial})
                    win = start_inds{trial}(sts):end_inds{trial}(sts);
                    extracted_sig = raw_sig_temp(win);
                    
                    raw_sig_temp(win) = raw_sig_temp(win) - avg_trial(1:length(win));
                    
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
        
end

fprintf(['-------Extracting data-------- \n \n'])

end