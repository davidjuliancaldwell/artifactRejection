function [processedSig,templateArray_cell,template] = template_subtract_iterative(raw_sig,varargin)
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

% default parameters
fs = 1.2207e04;
plotIt = 1;
pre = 0.4096;
post = 0.4096;
stimChans = [];

distanceDBscan = 0.001; % distance dbscan algorithm uses
distanceDBscan = 0.0005;
% define matrix of zeros
processedSig = zeros(size(raw_sig));
startInds = [];

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
        case 'startinds'
            startInds = varargin{i+1};
        case 'endinds'
            endInds = varargin{i+1};
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
lengthMax = 0; % length vector to build up the dictionary of templates later

% trim down start and end inds
startInds = cellfun(@(x) x+10,startInds,'un',0);
endInds = cellfun(@(x) x-15,endInds,'un',0);

%
for trial = 1:size(raw_sig,3)
    
    
    for chan = goodVec
        % get single trial
        raw_sig_temp = raw_sig(:,chan,trial);
        
        % default time to average
        default_time_average = endInds{trial}(1)-startInds{trial}(1)+1;
        
        
        
        lengthMax = max(default_time_average,lengthMax);
        
        avg_signal = zeros((endInds{trial}(1)-startInds{trial}(1)+1),length(startInds{trial}));
        for sts = 1:length(startInds{trial})
            win = startInds{trial}(sts):startInds{trial}(sts)+default_time_average-1;
            % avg_signal(:,sts) = raw_sig_temp(win);
            %avg_signal(:,sts) = raw_sig_temp(win) - raw_sig_temp(start_inds{trial}(sts));% take off first sample
            avg_signal(:,sts) = raw_sig_temp(win) - mean(raw_sig_temp(startInds{trial}(sts):startInds{trial}(sts)+presamps-3));% take off pre period
        end
        
        % find average stimulus
        avg_signal_mean = mean(avg_signal,2);
        
        if plotIt && (trial == 10 || trial == 1000)
            figure
            plot(raw_sig_temp,'linewidth',2)
            hold on
            plot(raw_sig(:,chan,trial),'linewidth',2)
            
            vline(startInds{trial})
            vline(endInds{trial},'g')
        end
        template{chan}{trial} = avg_signal;
    end
    
end
fprintf(['-------Finished getting artifacts-------- \n'])
%% build up dictionary

templateArray_cell = {};

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
    

    distanceDBscan = 0.75;
    %end
    
      %  distanceDBscan = 1e-4;
   % distanceDBscan = 7e-5;
    [c,ptsC,centres] = analyFunc.db_scan(templateArray,distanceDBscan,1,'cosine');
    
    vectorUniq = unique(ptsC);
    if plotIt
        figure
        hold on
        vectorUniq = unique(ptsC)';
        for i = 1:length(vectorUniq)
            plot(mean(templateArray(:,ptsC==i),2));
        end
    end
    templateArray_extracted = [];
    for i = 1:length(vectorUniq)
        meanTempArray = mean(templateArray(:,ptsC==i),2);
        % templateArray_extracted = [templateArray_extracted (meanTempArray )]; %no subtraction
          templateArray_extracted = [templateArray_extracted (meanTempArray-meanTempArray(1) )]; %take off first sample
        %templateArray_extracted = [templateArray_extracted (meanTempArray-median(meanTempArray(1:(presamps-6) )))]; % take off pre period
        
    end
    
    templateArray_cell{chan} = templateArray_extracted;
    
    fprintf(['-------Artifact Channel ' num2str(chan) ' -------- \n'])
    
    
end

fprintf(['-------Finished clustering artifacts-------- \n'])

%%
% now do the template subtraction
for trial = 1:size(raw_sig,3)
    
    for chan = goodVec
        
        raw_sig_temp = raw_sig(:,chan,trial);
        
        for sts = 1:length(startInds{trial})
            win = startInds{trial}(sts):startInds{trial}(sts)+default_time_average-1;
            extracted_sig = raw_sig_temp(win);
            
            % find best artifact
            templates = templateArray_cell{chan};
            correlations = corr(extracted_sig,templates);
            [~,index] = max((correlations));
            
            if max((correlations))>0.7
                template_subtract = templates(:,index);
            else
                template_subtract = zeros(size(templates(:,index)));
            end
            
            % scale
            %        span_art = max(template_subtract)-min(template_subtract);
            %         span_sig = max(extracted_sig)-min(extracted_sig);
            %         ratio = span_art/span_sig;
            %             ratio_max = max(template_subtract)/max(extracted_sig);
            %             ratio_min = min(template_subtract)/min(extracted_sig);
            %           ratio = mean([ratio_max,ratio_min]);
            %       template_subtract = template_subtract/ratio;
            
            raw_sig_temp(win) = raw_sig_temp(win) - template_subtract;
        end
        
        if plotIt && (trial == 10 || trial == 1000)
            figure
            plot(raw_sig_temp,'linewidth',2)
            hold on
            plot(raw_sig(:,chan,trial),'linewidth',2)
            
            vline(startInds{trial})
            vline(endInds{trial},'g')
        end
        
        processedSig(:,chan,trial) = raw_sig_temp;
        
    end
end
fprintf(['-------Extracting data-------- \n'])

end