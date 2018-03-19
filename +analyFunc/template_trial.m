function [processedSig,templateArrayCellOutput] = template_trial(templateTrial,rawSig,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;

addRequired(p,'template',@iscell);
addRequired(p,'rawSig',@isnumeric);

addParameter(p,'plotIt',0,@(x) x==0 || x ==1);

addParameter(p,'goodVec',[1:64],@isnumeric);
addParameter(p,'startInds',[],@iscell);
addParameter(p,'endInds',[],@iscell);


p.parse(templateTrial,rawSig,varargin{:});
templateTrial = p.Results.template;
plotIt = p.Results.plotIt;
rawSig = p.Results.rawSig;

goodVec = p.Results.goodVec;
startInds = p.Results.startInds;
endInds = p.Results.endInds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(['-------Trial Based Template-------- \n'])
templateArrayCellOutput = {};

for trial = 1:size(rawSig,3)
    
    for chan = goodVec
        rawSigTemp = rawSig(:,chan,trial);
        artifacts = templateTrial{chan}{trial};
        
        avgTrial = nanmean(artifacts,2);
        for sts = 1:length(startInds{trial}{chan})
            win = startInds{trial}{chan}(sts):endInds{trial}{chan}(sts);            
            rawSigTemp(win) = rawSigTemp(win) - avgTrial(1:length(win));
            
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