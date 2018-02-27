function [startInds,endInds] = get_artifact_indices(rawSig,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

validData = @(x) isnumeric(x) && size(x,3)>2;
addRequired(p,'rawSig',validData);

addParameter(p,'useFixedEnd',12207,@isnumeric);

addParameter(p,'pre',0.4096,@isnumeric);
addParameter(p,'plotIt',0,@(x) x==0 || x ==1);
addParameter(p,'post',0.4096,@isnumeric);
addParameter(p,'fixedDistance',2,@isnumeric);
addParameter(p,'fs',12207,@isnumeric);
addParameter(p,'goodVec',[1:64],@isnumeric);

p.parse(rawSig,varargin{:});

rawSig = p.Results.rawSig;
plotIt = p.Results.plotIt;
useFixedEnd = p.Results.useFixedEnd;
pre = p.Results.pre;
post = p.Results.post;
fixedDistance = p.Results.fixedDistance;
fs = p.Results.fs;
goodVec = p.Results.goodVec;


%% preprocess eco
presamps = round(pre/1e3 * fs); % pre time in sec

postsamps = round(post/1e3 * fs); %

fixedDistanceSamps = round(fixedDistance/1e3 * fs);

defaultWinAverage = fixedDistanceSamps ; %end_inds{trial}(1)-start_inds{trial}(1)+1;

% take diff of signal to find onset of stimulation train
diff_sig = permute(cat(3,zeros(size(rawSig,2), size(rawSig,3)),permute(diff(rawSig),[2 3 1])),[3 1 2]);

% find channel that has the max signal, and use this for subsequent
% analysis
[~,chanMax] = (max(max(diff_sig(:,goodVec,:))));
chanMax = chanMax(1);
lengthMax_vec = []; % length vector to build up the dictionary of templates later
fprintf(['-------Templates-------- \n'])

for trial = 1:size(rawSig,3)
    
    inds = find(abs(zscore(diff_sig(:,chanMax,trial)))>2);
    diff_bt_inds = [diff(inds)'];
    [~,inds_onset] = find(abs(zscore(diff_bt_inds))>2);
    startInds{trial} = [inds(1)-presamps; inds(inds_onset+1)-presamps];
    
    if useFixedEnd
        endInds{trial} = startInds{trial}+fixedDistanceSamps; % 17 to start
    else
        
        for idx = 1:length(startInds{trial})
            
            win = startInds{trial}(idx):startInds{trial}(idx)+defaultWinAverage; % get window that you know has the end of the stim pulse
            signal = rawSig(win,chanMax,trial);
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
            endInds{trial}(idx) = ct + startInds{trial}(idx) + postsamps;
            
        end
    end
    
    if plotIt
        figure
        plot(diff_sig(:,chanMax,trial))
        vline(startInds{trial})
        vline(endInds{trial},'g')
        
        figure
        plot(rawSig(:,chanMax,trial))
        vline(startInds{trial})
        vline(endInds{trial},'g')
    end
    
    
    
end
fprintf(['-------Finished getting artifacts-------- \n'])

end