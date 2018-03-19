function [processedSig,templateArrayCellOutput] = template_dictionary(templateArrayCell,templateTrial,rawSig,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;

addRequired(p,'templateArrayCell',@iscell);
addRequired(p,'templateTrial',@iscell);
addRequired(p,'rawSig',@isnumeric);

addParameter(p,'plotIt',0,@(x) x==0 || x ==1);

addParameter(p,'distanceMetricDbscan','eucl',@isstr);
addParameter(p,'distanceMetricSigMatch','eucl',@isstr);

addParameter(p,'goodVec',[1:64],@isnumeric);
addParameter(p,'startInds',[],@iscell);
addParameter(p,'endInds',[],@iscell);

addParameter(p,'recoverExp',1,@(x) x==0 || x ==1);


p.parse(templateArrayCell,templateTrial,rawSig,varargin{:});

templateArrayCell = p.Results.templateArrayCell;
templateTrial = p.Results.templateTrial;
rawSig = p.Results.rawSig;

plotIt = p.Results.plotIt;
distanceMetricDbscan = p.Results.distanceMetricDbscan;
distanceMetricSigMatch = p.Results.distanceMetricSigMatch;
goodVec = p.Results.goodVec;
startInds = p.Results.startInds;
endInds = p.Results.endInds;

recoverExp = p.Results.recoverExp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

templateArrayCellOutput = {};

fprintf(['-------Dictionary-------- \n'])

for chan = goodVec
    templateArrayExtracted = [];
    templateArray = templateArrayCell{chan};
    if strcmp(distanceMetricDbscan,'eucl')
        if max(templateArray(:)) < 3e-4
            distanceDBscan = 1e-3;
        else
            distanceDBscan = 1e-4;
        end
    else
        if max(templateArray(:)) < 3e-4
            distanceDBscan = 0.9;
        else
            distanceDBscan = 0.95;
        end
    end
    
    [c,ptsC,centres] = analyFunc.db_scan(templateArray,distanceDBscan,1,distanceMetricDbscan);
    
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
        templateArrayExtracted = [templateArrayExtracted (meanTempArray )]; %no subtraction
    end
    
    templateArrayCellOutput{chan} = templateArrayExtracted;
    
    fprintf(['-------Artifact Channel ' num2str(chan) ' -------- \n'])
    
    
end

fprintf(['-------Finished clustering artifacts-------- \n'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% now do the template subtraction
for trial = 1:size(rawSig,3)
    
    for chan = goodVec
        
        rawSigTemp = rawSig(:,chan,trial);
        templates = templateArrayCellOutput{chan};
        
        % add on the trial one
        templates = [templates mean(templateTrial{chan}{trial},2)];
        
        
        % ensure no subtraction of exponential
        
        if recoverExp
            templates = analyFunc.recover_EP(templates);
        end
        
        for sts = 1:length(startInds{trial}{chan})
            win = startInds{trial}{chan}(sts):endInds{trial}{chan}(sts);
            extractedSig = rawSigTemp(win);
            extractedSig = extractedSig - extractedSig(1);
            
            % find best artifact
            % get them to be the same length
            templatesSts = templates(1:length(extractedSig),:);
            
            switch distanceMetricSigMatch
                case 'correlation'
                    % correlation
                    correlations = corr(extractedSig,templatesSts);
                    [~,index] = max((correlations));
                    
                case 'cosineSim'
                    % - cosine similarity
                    denominator = sqrt(sum(extractedSig.*extractedSig)).*sqrt(sum(templatesSts.*templates));
                    numerator = (extractedSig'*templatesSts);
                    correlations = numerator./denominator;
                    [~,index] = max(abs(correlations));
                    
                case 'eucl'
                    % distance
                    v = templatesSts - repmat(extractedSig,1,size(templatesSts,2));
                    distance = sum(v.*v);
                    [~,index] = min(distance);
            end
            
            template_subtract = templatesSts(:,index);
            
            %             k = 2;
            %             all = templateArrayCell{chan};
            %             [u,s,v] = svd(all);
            %
            %             c = u(:,1:k);
            %             %c = template_subtract;
            %             a = rawSigTemp(win);
            %             d = (c'*c)\(c'*a);
            %             clean = a - (c'*c)\(c'*a);
            %
            rawSigTemp(win) = rawSigTemp(win) - template_subtract;
        end
        
        if plotIt && (trial == 10 || trial == 1000)
            figure
            plot(rawSigTemp,'linewidth',2)
            hold on
            plot(rawSig(:,chan,trial),'linewidth',2)
            
            vline(startInds{trial}{chan})
            vline(endInds{trial}{chan},'g')
            
            figure
            plot(extractedSig)
            hold on
            plot(template_subtract)
            plot(extractedSig-template_subtract)
            legend('extracted','template','subtracted');
        end
        
        processedSig(:,chan,trial) = rawSigTemp;
            fprintf(['-------Template Processed - Channel ' num2str(chan) '--' 'Trial ' num2str(trial) '-----\n'])

    end
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end