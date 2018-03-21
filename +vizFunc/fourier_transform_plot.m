function [f,P1] = fourier_transform_plot(fs_data,data,t)
% function to calculate and plot the results of applying the fourier
% transform to an epoched signal
% requires:
% fs_data - sampling rate in Hz
% data - time x channels signal of interest of interest
% t - a time vector for plotting (optional)

if ~exist('t','var')
    t = [0:size(data,1)-1];
end

[f,P1] = fourierTransformCalc(fs_data,data);

subplot(3,1,1)
plot((f),(P1),'linewidth',[2])
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([0 300])
ylim([0 2e-5])
set(gca,'fontsize',14)
hold on

subplot(3,1,2)
loglog((f),(P1),'linewidth',[2])
title('Single-Sided Amplitude Spectrum of X(t) -loglog')
xlabel('f (Hz)')
ylabel('|P1(f)|')
% xlim([0 500])
ylim([10e-10 10e-4])

set(gca,'fontsize',14)
hold on

subplot(3,1,3)
plot(t,data,'linewidth',[2])
title('time series')
xlabel('time (ms)')
ylabel('ampltitude (V)')
set(gca,'fontsize',14)
xlim([-500 500])
ylim ([-1e-4 1e-4])
hold on


end