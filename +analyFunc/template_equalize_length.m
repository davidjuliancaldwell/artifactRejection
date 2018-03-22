function [template,templateArrayCell] = template_equalize_length(templateCell,rawSig,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;

addRequired(p,'templateCell',@iscell);
addRequired(p,'rawSig',@isnumeric);

addParameter(p,'goodVec',[1:64],@isnumeric);
addParameter(p,'startInds',[],@iscell);
addParameter(p,'lengthMax',25,@isnumeric);

p.parse(templateCell,rawSig,varargin{:});
templateCell = p.Results.templateCell;
rawSig = p.Results.rawSig;
goodVec = p.Results.goodVec;
startInds = p.Results.startInds;
lengthMax = p.Results.lengthMax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for chan = goodVec
    templateArray = [];
    templateArrayExtracted = [];
    lengthMaxChan = lengthMax(chan);
    
    for trial = 1:size(rawSig,3)
        artifactsCell = templateCell{chan}{trial};
        artifactsMat = [];
        for sts = 1:length(startInds{trial}{chan})
            
            artifactsTrial = artifactsCell{sts};
            
            if size(artifactsTrial,1) < lengthMaxChan
                amntPad = lengthMaxChan - size(artifactsTrial,1);
                artifacts_pad = padarray(artifactsTrial,amntPad,0,'post');
            else
                artifacts_pad = artifactsTrial;
            end
            
            artifactsMat(:,sts) = artifacts_pad;
            
        end
        template{chan}{trial} = artifactsMat;
        templateArray = [templateArray artifactsMat];
        templateArrayCell{chan} = templateArray;
    end
end

fprintf(['-------Finished making artifacts the same length-------- \n'])

end