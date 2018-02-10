function [ wavelet_norm ] = normalize_spectrogram(data,t)
%NORMALIZE_SPECTROGRAM Summary of this function goes here
%   Detailed explanation goes here

pre_data = data(:,t<0,:,:);
post_data = data(:,t>0,:,:);
wavelet_norm = zeros(size(data));


for channel = 1:size(data,3)
    for trial = 1:size(data,4)
        wavelet_norm(:,:,channel,trial) = normalize_data(data(:,:,channel,trial),pre_data(:,:,channel,trial));
    end
end


end

