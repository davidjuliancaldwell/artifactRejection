function [normalizedData] = normalize_spectrogram_wavelet_on_avg_signal(dataRef,data)

% freq x time x channel  

dataRefAvg = squeeze(mean(dataRef,4));
dataAvg = squeeze(mean(data,4));

normalizedData = zeros(size(dataAvg));

for chan = 1:size(data,3)
    
        refTempTrial = squeeze(dataRefAvg(:,:,chan));
        mRef1=squeeze(mean(refTempTrial,2));
        sRef1=squeeze(std(refTempTrial'))';
        
        dataTempTrial = squeeze(dataAvg(:,:,chan));
        normalizedData(:,:,chan)=(dataTempTrial-mRef1*ones(1,size(dataTempTrial,2)))./(sRef1*ones(1,size(dataTempTrial,2)));
    
end

end

