function [normalizedData] = normalize_spectrogram_wavelet_fused_baseline(dataRef,data)

% time x freq x channel x trial

normalizedData = zeros(size(data));

for chan = 1:size(data,3)
    
    trial_baseline = squeeze(dataRef(:,:,chan,:));
    size_trial_baseline = size(trial_baseline);
    
    trial_baseline = reshape(trial_baseline,size_trial_baseline(1),[],1);
        
    mRef1=squeeze(mean(trial_baseline,2));
    sRef1=squeeze(std(trial_baseline'))';
    
    
    for trial = 1:size(data,4)
        
        dataTempTrial = squeeze(data(:,:,chan,trial));
        normalizedData(:,:,chan,trial)=(dataTempTrial-mRef1*ones(1,size(dataTempTrial,2)))./(sRef1*ones(1,size(dataTempTrial,2)));
    end
end

end

