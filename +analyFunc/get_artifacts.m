function [lengthMax,templateCell] = get_artifacts(rawSig,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

validData = @(x) isnumeric(x) && size(x,3)>2;
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

lengthMaxVec = [];

for trial = 1:size(rawSig,3)
    
    for chan = goodVec
        % get single trial
        rawSigTemp = rawSig(:,chan,trial);
        
        avgSignal = {};
        lengthMaxVecTemp = [];
        for sts = 1:length(startInds{trial}{chan})
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
        
        lengthMaxVec = [lengthMaxVec max(lengthMaxVecTemp)];
        templateCell{chan}{trial} = avgSignal;
        
        if plotIt && (trial == 10 || trial == 1000)
            figure
            plot(rawSigTemp,'linewidth',2)
            vline(startInds{trial})
            vline(endInds{trial},'g')
        end
        
    end
    
end

lengthMax = max(lengthMaxVec);

end