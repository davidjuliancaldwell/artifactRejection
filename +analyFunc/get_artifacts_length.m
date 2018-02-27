function [lengthMaxVec,templateCell] = get_artifacts(rawSig,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

validData = @(x) isnumeric(x) && size(x,3)>2;
addRequired(p,'rawSig',validData);
addParameter(p,'plotIt',0,@(x) x==0 || x ==1);
addParameter(p,'goodVec',[1:64],@isnumeric);
addParameter(p,'startInds',[],@isnumeric);
addParameter(p,'endInds',[],@isnumeric);



p.parse(rawSig,varargin{:});
rawSig = p.Results.rawSig;
startInds = p.Results.startInds;
endInds = p.Results.endInds;
goodVec = p.Results.goodVec;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for trial = 1:size(rawSig,3)
    
    for chan = goodVec
        % get single trial
        rawSigTemp = rawSig(:,chan,trial);
        
        avgSignal = {};
        lengthMaxVecTemp = [];
        for sts = 1:length(startInds{trial})
            win = startInds{trial}(sts):endInds{trial}(sts);
            lengthMax_temp = length(win);
            avgSignal{sts} = rawSigTemp(win) - mean(rawSigTemp(startInds{trial}(sts):startInds{trial}(sts)+presamps-8));% take off average of first 3 samples
            
            lengthMaxVecTemp = [lengthMaxVecTemp lengthMax_temp];
            
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
    lengthMax = max(lengthMaxVec);

end

end