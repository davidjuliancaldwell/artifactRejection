%% script to create parts of figure 2 - need processed sig


%%

avgResponse = mean(processedSig,3);
avgRaw = mean(dataInt,3);

% look at the FFT difference
[f,P1] = helpFunc.fourier_transform_calc(fsData,avgResponse);
[fRaw,P1Raw] = helpFunc.fourier_transform_calc(fsData,avgRaw);

idx = 28;
figure
plot(fRaw,P1Raw(:,idx),'Color',[0.5 0.5 0.5],'LineWidth',2)
ylabel('Log Power')
xlabel('Frequency (Hz)')
set(gca,'fontsize',18)
title('Frequency Spectrum of Raw Signal')
xlim([0 700])
figure

semilogy(fRaw,P1Raw(:,idx),'Color',[0.5 0.5 0.5],'LineWidth',2)
ylabel('Log Power')
xlabel('Frequency (Hz)')
set(gca,'fontsize',18)
title('Frequency Spectrum of Raw Signal')
xlim([0 700])

idx = 28;
figure
plot(f,P1(:,idx),'Color',[0.5 0.5 0.5],'LineWidth',2)
ylabel('Log Power')
xlabel('Frequency (Hz)')
set(gca,'fontsize',18)
title('Frequency Spectrum of Processed Signal')
xlim([0 700])

figure

semilogy(f,P1(:,idx),'Color',[0.5 0.5 0.5],'LineWidth',2)
ylabel('Log Power')
xlabel('Frequency (Hz)')
set(gca,'fontsize',18)
title('Frequency Spectrum of Processed Signal')

xlim([0 700])
