function multiple_visualizations(processedSig,rawSig,varargin)
% function to run multiple visualizations and analyses to assess the quality
% of the processing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get inputs
p = inputParser;

validData = @(x) isnumeric(x) && size(x,3)>2;
addRequired(p,'processedSig',validData);
addRequired(p,'rawSig',validData);

addParameter(p,'type','linear',@isstr);

addParameter(p,'xlims',[-100 1000],@isnumeric);
addParameter(p,'ylims',[-500 500],@isnumeric);
addParameter(p,'trainDuration',[0 400],@isnumeric);
addParameter(p,'tEpoch',0.2,@isnumeric);
addParameter(p,'stimChans',[],@isnumeric);
addParameter(p,'bads',[],@isnumeric);
addParameter(p,'chanIntList',[1,2,3],@isnumeric);
addParameter(p,'fs',12207,@isnumeric);
addParameter(p,'templateTrial',@iscell)
addParameter(p,'templateDictCell',@iscell);
addParameter(p,'modePlot','avg',@isstr);

p.parse(processedSig,rawSig,varargin{:});

rawSig = p.Results.rawSig;

type = p.Results.type;

stimChans = p.Results.stimChans;
bads = p.Results.bads;
fs = p.Results.fs;

templateTrial = p.Results.templateTrial;
templateDictCell = p.Results.templateDictCell;
trainDuration = p.Results.trainDuration;
xlims = p.Results.xlims;
ylims = p.Results.ylims;
tEpoch = p.Results.tEpoch;
chanIntList = p.Results.chanIntList;
modePlot = p.Results.modePlot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
numChans = size(rawSig,2);
[goods,goodVec] = helpFunc.goodChannel_extract('bads',bads,'stimchans',stimChans,'numChans',numChans);

if (strcmp(type,'dictionary') || strcmp(type,'trial') || strcmp(type,'average')) && (exist('templateDictCell','var') || exist('templateTrial','var'))

    if exist('templateTrial','var')
    figure
    for j = goodVec
        subplot(8,8,j)
        hold on
        for i = 1:size(templateTrial{j},2)
            timeVec = [0:size(templateTrial{j}{:,i},1)-1];
            
            vizFunc.plotBTLError(timeVec,templateTrial{j}{:,i},'CI',rand(1,3));
        end
        title(['Channel ' num2str(j)])
    end
    
    % plot the average template dictionary if using a dictionary method
    if strcmp(type,'dictionary') && exist('templateDictCell','var')
        figure
        hold on
        for j = goodVec
            subplot(8,8,j)
            hold on
            for i = 1:size(templateDictCell{j},2)
                timeVec = [0:size(templateDictCell{j},1)-1];
                
                plot(timeVec,templateDictCell{j}(:,i),'linewidth',2);
            end
            title(['Channel ' num2str(j)])
        end
    end
end

avgResponse = mean(processedSig,3);
avgRaw = mean(rawSig,3);

vizFunc.smallMultiples_artifactReject(processedSig,tEpoch,'type1',stimChans,'type2',0,'modePlot',modePlot,'highlightRange',trainDuration)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make a list of the channels of interest to visualize one at a time

for ind = chanIntList
    
    exampChan = mean(squeeze(processedSig(:,ind,:)),2);
    
    figure
    ax1 = subplot(2,1,1);
    plot(1e3*tEpoch,1e6*exampChan,'linewidth',2);
    xlim(xlims)
    ylim(ylims)
    
    title(['Processed Signal - Channel ' num2str(ind)])
    clear exampChan
    
    
    ax2 = subplot(2,1,2);
    exampChan = mean(squeeze(rawSig(:,ind,:)),2);
    plot(1e3*tEpoch,1e6*exampChan,'linewidth',2);
    xlim(xlims)
    ylim(ylims)
    xlabel('time (ms)')
    ylabel('Voltage (\muV)')
    title(['Raw Signal Average - Channel ' num2str(ind)])
    
    linkaxes([ax1,ax2],'xy')
    
    clear exampChan
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% look at the FFT difference
[f,P1] = helpFunc.fourierTransformCalc(fs,avgResponse);
[fRaw,P1Raw] = helpFunc.fourierTransformCalc(fs,avgRaw);

vizFunc.smallMultiples_fourier(P1Raw,fRaw,'type1',stimChans,'type2',0)
legend('raw')
vizFunc.smallMultiples_fourier(P1,f,'type1',stimChans,'type2',0,'newfig',0)
legend('processed')

vizFunc.smallMultiples_fourier(P1Raw,fRaw,'type1',stimChans,'type2',0,'plotLog',1)
legend('raw')
vizFunc.smallMultiples_fourier(P1,f,'type1',stimChans,'type2',0,'newfig',0,'plotLog',1)
legend('processed')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% look at the RMS reduction within particular section
processedSigRms = helpFunc.rms_func(avgResponse(tEpoch<1000 & tEpoch>-100,:));
rawSigRms = helpFunc.rms_func(avgRaw(tEpoch<1000 & tEpoch>-100,:));

% look at decibel reduction for each channel

rms_db = 20*log10(processedSigRms./rawSigRms);

figure
numBins = 10;
histogram(rms_db,numBins)
title('RMS Reduction in Decibels')
xlabel('magnitude of decibel decrease')
ylabel('count')
set(gca,'fontsize',14)

end