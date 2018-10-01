function [templateCell,lengthMax,maxAmps,maxLocation] = get_artifacts(rawSig,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

validData = @(x) isnumeric(x);
addRequired(p,'rawSig',validData);
addParameter(p,'plotIt',0,@(x) x==0 || x ==1);
addParameter(p,'goodVec',[1:64],@isnumeric);
addParameter(p,'startInds',[],@iscell);
addParameter(p,'endInds',[],@iscell);
addParameter(p,'normalize','firstSamp',@isstr);
addParameter(p,'amntPreAverage',12,@isnumeric)

p.parse(rawSig,varargin{:});
rawSig = p.Results.rawSig;
plotIt = p.Results.plotIt;
startInds = p.Results.startInds;
endInds = p.Results.endInds;
goodVec = p.Results.goodVec;
normalize = p.Results.normalize;
amntPreAverage = p.Results.amntPreAverage;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lengthMaxVecTrial = zeros(size(rawSig,3),size(rawSig,2));

for trial = 1:size(rawSig,3) % loop through trials 
    lengthMaxChan = zeros(1,size(rawSig,2));
    
    for chan = goodVec % loop through channels
        % get single trial
        rawSigTemp = rawSig(:,chan,trial);
        
        avgSignal = {};
        lengthMaxVecTemp = [];
        for sts = 1:length(startInds{trial}{chan}) % loop through individual stimulation epochs 
            win = startInds{trial}{chan}(sts):endInds{trial}{chan}(sts);
            lengthMaxTemp = length(win);
            
            switch normalize
                case 'preAverage'
                    avgSignal{sts} = rawSigTemp(win) - mean(rawSigTemp(startInds{trial}{chan}(sts):startInds{trial}{chan}(sts)+amntPreAverage));% take off average of first x samples
                case 'none'
                    avgSignal{sts} = rawSigTemp(win);
                case 'firstSamp'
                    avgSignal{sts} = rawSigTemp(win) - rawSigTemp(startInds{trial}{chan}(sts));
            end
            
            lengthMaxVecTemp = [lengthMaxVecTemp lengthMaxTemp];
            
        end
        
        lengthMaxChan(chan) = max(lengthMaxVecTemp);
        
        templateCell{chan}{trial} = avgSignal;
        
        if plotIt && (trial == 10 || trial == 1000)
            figure
            plot(rawSigTemp,'linewidth',2)
            vline(startInds{trial})
            vline(endInds{trial},'g')
        end
        
    end
    lengthMaxVecTrial(trial,:) = lengthMaxChan; 
    
end

% figure out the maximum amplitude artifact for each given channel and
% trial 

maxAmps = squeeze(max(rawSig,[],1));

% find maximum index for reducing dimensionality later 
[~,maxChan] = max(maxAmps,[],1);
maxChan = maxChan(1);
[~,maxTrial] = max(maxAmps,[],2);
maxTrial = maxTrial(maxChan);
[~,maxLocation] = cellfun(@max,(templateCell{maxChan}{maxTrial}));
maxLocation = median(maxLocation);

% get the maximum length of any given artifact window for each channel
lengthMax = max(lengthMaxVecTrial,[],1);

end