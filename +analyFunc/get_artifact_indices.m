function [startInds,endInds] = get_artifact_indices(rawSig,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

validData = @(x) isnumeric(x) && size(x,3)>2;
addRequired(p,'rawSig',validData);

addParameter(p,'useFixedEnd',0,@(x) x==0 || x ==1);

addParameter(p,'pre',0.4096,@isnumeric);
addParameter(p,'plotIt',0,@(x) x==0 || x ==1);
addParameter(p,'post',0.4096,@isnumeric);
addParameter(p,'fixedDistance',2,@isnumeric);
addParameter(p,'fs',12207,@isnumeric);
addParameter(p,'goodVec',[1:64],@isnumeric);
addParameter(p,'chanInt',15,@isnumeric);

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


%% preprocess eco
presamps = round(pre/1e3 * fs); % pre time in sec

postsamps = round(post/1e3 * fs); %

fixedDistanceSamps = round(fixedDistance/1e3 * fs);

defaultWinAverage = fixedDistanceSamps ; %end_inds{trial}(1)-start_inds{trial}(1)+1;

% take diff of signal to find onset of stimulation train
diffSig = permute(cat(3,zeros(size(rawSig,2), size(rawSig,3)),permute(diff(rawSig),[2 3 1])),[3 1 2]);

% find channel that has the max signal, and use this for subsequent
% analysis
[~,chanMax] = (max(max(diffSig(:,goodVec,:))));
chanMax = chanMax(1);
lengthMax_vec = []; % length vector to build up the dictionary of templates later
fprintf(['-------Templates-------- \n'])

if ~exist('chanInt','var')
    chanInt = chanMax;
end

for trial = 1:size(rawSig,3)
    
    inds = find(abs(zscore(diffSig(:,chanMax,trial)))>2);
    diffBtInds = [diff(inds)'];
    [~,indsOnset] = find(abs(zscore(diffBtInds))>2);
    
    for chan = goodVec
        
        startInds{trial}{chan} = [inds(1)-presamps; inds(indsOnset+1)-presamps];
        
        if useFixedEnd
            endInds{trial}{chan} = startInds{trial}{chan}+fixedDistanceSamps; % 17 to start
        else
            
            for idx = 1:length(startInds{trial}{chan})
                
                win = startInds{trial}{chan}(idx):startInds{trial}{chan}(idx)+defaultWinAverage; % get window that you know has the end of the stim pulse
                signal = rawSig(win,chan,trial);
                diff_signal = diffSig(win,chan,trial);
                
                last = find(abs(zscore(signal))>0.2,1,'last');
                last2 = find(abs(zscore(diff_signal))>5e-3,1,'last')+1;
                ct = max(last, last2);
                
                if length(win) - ct > 8
                    
                    x = [ct:length(win)]';
                    y = signal(x);
                    [f2,gof,output] = fit(x,y,'exp2');
                    func_fit = @(x) f2.a*exp(f2.b*x) + f2.c*exp(f2.d*x);
                    
                    % if its a good fit, set the end index to include that
                    
                    if gof.adjrsquare>0.8
                        ct = length(win);
                    end
                    
                end
                
                endInds{trial}{chan}(idx) = ct + startInds{trial}{chan}(idx) + postsamps;
                
            end
        end
        
    end
    
    
    if plotIt
        figure
        plot(diffSig(:,chanInt,trial))
        vline(startInds{trial}{chanInt})
        vline(endInds{trial}{chanInt},'g')
        
        figure
        plot(rawSig(:,chanInt,trial))
        vline(startInds{trial}{chanInt})
        vline(endInds{trial}{chanInt},'g')
    end
    
    fprintf(['-------Finished getting artifacts - Trial ' num2str(trial) '-------- \n'])
    
end