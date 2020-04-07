function [normalizedData] = normalize_spectrogram_wavelet(dataRef,data)

% freq x time x channel x trial 

normalizedData = zeros(size(data));

for chan = 1:size(data,3)

    
    for trial = 1:size(data,4)
        refTempTrial = squeeze(dataRef(:,:,chan,trial));
        mRef1=squeeze(mean(refTempTrial,2));
        sRef1=squeeze(std(refTempTrial'))';
        
        dataTempTrial = squeeze(data(:,:,chan,trial));
        normalizedData(:,:,chan,trial)=(dataTempTrial-mRef1*ones(1,size(dataTempTrial,2)))./(sRef1*ones(1,size(dataTempTrial,2)));
    end
end

end

