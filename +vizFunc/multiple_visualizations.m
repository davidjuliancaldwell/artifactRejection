function multiple_visualizations(processedSig,rawSig,varargin)
% function to run multiple visualizations and analyses to assess the quality
% of the processing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check some basic data requirements
if nargin == 0
    error ('You must supply data');
end

if length (size (rawSig)) > 3
    error ('Input data can not have more than three dimensions.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% defaults
type = 'linear';
template = {};
templateDictCell = {};
chanInt = [];
chanIntList = [];
stimChans = [];
t_epoch = [0:size(rawSig,1)-1];
xlims = [-100 1000];
ylims = [-500 500];
trainDuration = [0 400];
fs_data = 12207;

% variable input
for i=1:2:(length(varargin)-1)
    if ~ischar (varargin{i})
        error (['Unknown type of optional parameter name (parameter' ...
            ' names must be strings).']);
    end
    % change the value of parameter
    switch lower (varargin{i})
        case 'type'
            type= varargin{i+1};
        case 'chanint'
            chanInt = varargin{i+1};
        case 'chanintlist'
            chanIntList = varargin{i+1};
        case 'stimchans'
            stimChans = varargin{i+1};
        case 't_epoch'
            t_epoch = varargin{i+1};
        case 'trainduration'
            trainDuration = varargin{i+1};
        case 'xlims'
            xlims = varargin{i+1};
        case 'ylims'
            ylims = varargin{i+1};
        case 'fs_data'
            fs_data = varargin{i+1};
        case 'template'
            template = varargin{i+1};
        case 'templatedictcell'
            templateDictCell = varargin{i+1};
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
numChans = size(rawSig,2);
[goods,goodVec] = helpFunc.goodChannel_extract('stimchans',stimChans,'numChans',numChans);

if strcmp(type,'dictionary') || strcmp(type,'trial') || strcmp(type,'average')

    figure
    for j = goodVec
        subplot(8,8,j)
        hold on
        for i = 1:size(template{j},2)
            timeVec = [0:size(template{j}{:,i},1)-1];
            
            vizFunc.plotBTLError(timeVec,template{j}{:,i},'CI',rand(1,3)');
        end
        title(['Channel ' num2str(j)])
    end
    
    % plot the average template dictionary if using a dictionary method
    if strcmp(type,'dictionary')
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

vizFunc.smallMultiples_artifactReject(avgResponse,t_epoch,'type1',stimChans,'type2',0,'average',1,'highlight_range',trainDuration)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make a list of the channels of interest to visualize one at a time

for ind = chanIntList
    
    exampChan = mean(squeeze(processedSig(:,ind,:)),2);
    
    figure
    ax1 = subplot(2,1,1);
    plot(1e3*t_epoch,1e6*exampChan,'linewidth',2);
    xlim(xlims)
    ylim(ylims)
    
    title(['Processed Signal - Channel ' num2str(ind)])
    clear exampChan
    
    
    ax2 = subplot(2,1,2);
    exampChan = mean(squeeze(rawSig(:,ind,:)),2);
    plot(1e3*t_epoch,1e6*exampChan,'linewidth',2);
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
[f,P1] = helpFunc.fourierTransformCalc(fs_data,avgResponse);
[fRaw,P1Raw] = helpFunc.fourierTransformCalc(fs_data,avgRaw);

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
processedSig_rms = helpFunc.rms_func(avgResponse(t_epoch<1000 & t_epoch>-100,:));
rawSig_rms = helpFunc.rms_func(avgRaw(t_epoch<1000 & t_epoch>-100,:));

% look at decibel reduction for each channel

rms_db = 20*log10(processedSig_rms./rawSig_rms);

figure
numBins = 10;
histogram(rms_db,numBins)
title('RMS Reduction in Decibels')
xlabel('magnitude of decibel decrease')
ylabel('count')
set(gca,'fontsize',14)

end