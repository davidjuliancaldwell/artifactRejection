function [processedSig,templateArrayCellOutput] = template_average(templateArrayCell,rawSig,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;

addRequired(p,'templateArrayCell',@iscell);
addRequired(p,'rawSig',@isnumeric);

addParameter(p,'plotIt',0,@(x) x==0 || x ==1);

addParameter(p,'goodVec',[1:64],@isnumeric);
addParameter(p,'startInds',[],@iscell);
addParameter(p,'endInds',[],@iscell);


p.parse(templateArrayCell,rawSig,varargin{:});
templateArrayCell = p.Results.templateArrayCell;
plotIt = p.Results.plotIt;
rawSig = p.Results.rawSig;

goodVec = p.Results.goodVec;
startInds = p.Results.startInds;
endInds = p.Results.endInds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

templateArrayCellOutput = {};
processedSig = zeros(size(rawSig));

fprintf(['-------Average Template-------- \n'])

for chan = goodVec
    templateArray = templateArrayCell{chan};
    
    avgSignalMean = nanmean(templateArray,2);
    templateArrayCell{chan} = avgSignalMean;
    
    for trial = 1:size(rawSig,3)
        rawSigTemp = rawSig(:,chan,trial);
        
        for sts = 1:length(startInds{trial})
            win = startInds{trial}{chan}(sts):endInds{trial}{chan}(sts);
            rawSigTemp(win) = rawSigTemp(win) - avgSignalMean(1:length(win));
        end
        
        if plotIt && (trial == 10 || trial == 1000)
            figure
            plot(rawSigTemp,'linewidth',2)
            hold on
            plot(rawSig(:,chan,trial),'linewidth',2)
            
            vline(startInds{trial}{chan})
            vline(endInds{trial}{chan},'g')
            
        end
        
        processedSig(:,chan,trial) = rawSigTemp;
    end
    
end

end