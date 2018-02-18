function [processedSig,templateArray_cell,template,start_inds,end_inds] = templateSubtract(raw_sig,varargin)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

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

for trial = 1:size(raw_sig,3)
    
    % find onset
    
    %             foo = mean(temp,2);
    %         lastsample = round(0.040 * efs);
    %         foo(lastsample:end) = foo(lastsample-1);
    %
    %         last = find(abs(zscore(foo))>1,1,'last');
    %         last2 = find(abs(diff(foo))>30e-6,1,'last')+1;
    %
    %         zc = false;
    %
    %         if (isempty(last2))
    %             if (isempty(last))
    %                 error ('something seems wrong in the triggered average');
    %             else
    %                 ct = last;
    %             end
    %         else
    %             if (isempty(last))
    %                 ct = last2;
    %             else
    %                 ct = max(last, last2);
    %             end
    %         end
    %         % try getting rid of this part for 0b5a2e to conserve that initial
    %         % spike DJC 1-7-2016
    %         while (~zc && ct <= length(foo))
    %             zc = sign(foo(ct-1)) ~= sign(foo(ct));
    %             ct = ct + 1;
    %         end
    %         % consider 3 ms? DJC - 1-5-2016
    %         if (ct > max(last, last2) + 0.10 * efs) % marched along more than 10 msec, probably gone to far
    %             ct = max(last, last2);
    %         end
    %
    %         %
    
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
        
        avg_signal = zeros((end_inds{trial}(1)-start_inds{trial}(1)+1),length(start_inds{trial}));
        for sts = 1:length(start_inds{trial})
            win = start_inds{trial}(sts):start_inds{trial}(sts)+default_time_average-1;
            % avg_signal(:,sts) = raw_sig_temp(win);
            %avg_signal(:,sts) = raw_sig_temp(win) - raw_sig_temp(start_inds{trial}(sts));% take off first sample
            avg_signal(:,sts) = raw_sig_temp(win) - mean(raw_sig_temp(start_inds{trial}(sts):start_inds{trial}(sts)+presamps-8));% take off pre period
        end
        
        % find average stimulus
        avg_signal_mean = mean(avg_signal,2);
        
        if plotIt && (trial == 10 || trial == 1000)
            figure
            plot(raw_sig_temp,'linewidth',2)
            hold on
            plot(raw_sig(:,chan,trial),'linewidth',2)
            
            vline(start_inds{trial})
            vline(end_inds{trial},'g')
        end
        template{chan}{trial} = avg_signal;
    end
    
end
fprintf(['-------Finished getting artifacts-------- \n'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% build up dictionary

templateArray_cell = {};

switch type
    case 'dictionary'
        for chan = goodVec
            templateArray = [];
            templateArray_extracted = [];
            
            for trial = 1:size(raw_sig,3)
                artifacts = template{chan}{trial};
                if size(artifacts,1) < lengthMax
                    amntPad = lengthMax - size(artifacts,1);
                    artifacts = padarray(artifacts,[amntPad 0],'post');
                end
                templateArray = [templateArray artifacts];
            end
            
            % adaptively change dbscan parameter
            %     if max(templateArray(:)) < 3e-4
            %         distanceDBscan = 0.0001;
            %     else
            %         distanceDBscan = 0.001;
            %     end
            
            %distanceDBscan = 0;
            %  if max(templateArray(:)) < 3e-4
            % distanceDBscan = 0.99;
            %  else
            distanceDBscan = 0.95;
            % end
            %  if chan == 12
            %         distanceDBscan = 0.995;
            % else
            % distanceDBscan = 0.99;
            
            %end
            distanceDBscan = 1e-7;
            distanceDBscan = 0.95;
            
            if max(templateArray(:)) < 3e-4
                distanceDBscan = 0.9;
            else
                distanceDBscan = 0.95;
            end
            
            %[c,ptsC,centres] = dbscan(templateArray,distanceDBscan,1,'corr');
            [c,ptsC,centres] = dbscan(templateArray,distanceDBscan,1,'corr');
            
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
                %  templateArray_extracted = [templateArray_extracted (meanTempArray-meanTempArray(1) )]; %take off first sample
                % templateArray_extracted = [templateArray_extracted (meanTempArray-median(meanTempArray(1:(presamps-6) )))]; % take off pre period
                
            end
            
            templateArray_cell{chan} = templateArray_extracted;
            
            fprintf(['-------Artifact Channel ' num2str(chan) ' -------- \n'])
            
            
        end
        
        fprintf(['-------Finished clustering artifacts-------- \n'])
        
        % now do the template subtraction
        for trial = 1:size(raw_sig,3)
            
            for chan = goodVec
                
                raw_sig_temp = raw_sig(:,chan,trial);
                
                for sts = 1:length(start_inds{trial})
                    win = start_inds{trial}(sts):start_inds{trial}(sts)+default_time_average-1;
                    extracted_sig = raw_sig_temp(win);
                    
                    % find best artifact
                    
                    templates = templateArray_cell{chan};
                    
                    % add on the trial one
                    templates = [templates mean(template{chan}{trial},2)];
                    
                    % correlation
                    
                    correlations = corr(extracted_sig,templates);
                    [~,index] = max((correlations));
                    template_subtract = templates(:,index);
                    
                    % - cosine similarity
                    %             denominator = sqrt(sum(extracted_sig.*extracted_sig)).*sqrt(sum(templates.*templates));
                    %             numerator = (extracted_sig'*templates);
                    %             correlations = numerator./denominator;
                    %                         [~,index] = max(abs(correlations));
                    
                    % distance
                    
                    v = templates - repmat(extracted_sig,1,size(templates,2));
                    distance = sum(v.*v);
                    [~,index] = min(distance);
                    
                    template_subtract = templates(:,index);
                    
                    %             if max((correlations))>0.9
                    %                 template_subtract = templates(:,index);
                    %             else
                    %                 template_subtract = zeros(size(templates(:,index)));
                    %             end
                    %
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
                    
                    vline(start_inds{trial})
                    vline(end_inds{trial},'g')
                    
                    figure
                    plot(extracted_sig)
                    hold on
                    plot(template_subtract)
                    plot(extracted_sig-template_subtract)
                    legend('extracted','template','subtracted');
                    %
                    %
                    %             plot(a)
                    %             plot(extracted_sig-a)
                    %             legend('extracted','template','subtracted','trial_temp','subt_trial');
                    
                end
                
                processedSig(:,chan,trial) = raw_sig_temp;
                
            end
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    case 'trial'
        for trial = 1:size(raw_sig,3)
            
            for chan = goodVec
                templateArray = [];
                raw_sig_temp = raw_sig(:,chan,trial);
                artifacts = template{chan}{trial};
                
                
                foo = mean(artifacts,2);
                %         last = find(abs(zscore(foo))>1,1,'last');
                %         last2 = find(abs(diff(foo))>20e-6,1,'last')+1;
                %         last3 = find(abs(diff(foo))>20e-6,1,'last')+1;
                %
                %         x = [last2:length(foo)]';
                %         y = foo(x);
                %         [f2,gof,output] = fit(x,y,'exp2');
                %         func_fit = @(x) f2.a*exp(f2.b*x) + f2.c*exp(f2.d*x);
                %
                %         if gof.adjrsquare>0.9
                %             foo(x) = func_fit(x);
                %         end
                %
                %         if plotIt
                %             figure
                %             plot(foo)
                %             vline(last2)
                %             vline(last3,'g')
                %
                %             figure
                %             plot(foo)
                %             hold on
                %             plot(x,y,'linewidth',2)
                %         end
                %
                for sts = 1:length(start_inds{trial})
                    win = start_inds{trial}(sts):start_inds{trial}(sts)+default_time_average-1;
                    extracted_sig = raw_sig_temp(win);
                    
                    raw_sig_temp(win) = raw_sig_temp(win) - foo;
                    
                end
                
                processedSig(:,chan,trial) = raw_sig_temp;         
            end
        end
             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    case 'average'
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
        templateArray_cell = {}
end

fprintf(['-------Extracting data-------- \n'])

end