function [f,P1] = fourier_transform_calc(fs_data,data)
% function to compute the fourier transform of a time x channels signal
% requires:
% fs_data - sampling frequency in Hz
% data - samples x channels 

fs = fs_data;
T = 1/fs;
L = size(data,1);

Y = fft(data);
P2 = abs(Y/L);
P1 = P2(1:floor(L/2)+1,:);
P1(2:end-1,:) = 2*P1(2:end-1,:);
f = fs*(0:floor(L/2))/L;

end