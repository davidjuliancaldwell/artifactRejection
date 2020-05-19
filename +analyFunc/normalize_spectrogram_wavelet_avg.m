function [normalizedData] = normalize_spectrogram_wavelet_avg(dataRef,data)

% time x freq x channel x trial

normalizedData = zeros(size(data));

for chan = 1:size(data,3)
    
    trial_mean = mean(squeeze(dataRef(:,:,chan,:)),3);
        
    mRef1=squeeze(mean(trial_mean,2));
    sRef1=squeeze(std(trial_mean'))';
    
    
    for trial = 1:size(data,4)
        
        dataTempTrial = squeeze(data(:,:,chan,trial));
        normalizedData(:,:,chan,trial)=(dataTempTrial-mRef1*ones(1,size(dataTempTrial,2)))./(sRef1*ones(1,size(dataTempTrial,2)));
    end
end

end

