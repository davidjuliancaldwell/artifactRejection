% inspired by Trebaul et al. 2016
% RC circuit

% this is the sampling rate to make the artifacts at!
%fs = 12207*5;
fs = 100000; % change this as desired!!!

% this is the duration of the pulse
durPulse = 0.0012; % change this as desired!!!

% this is the desired frequency to downsample to
fDesired = 12207; % change this as desired!!!

% this is the vector of samples for shifting
circshiftIndices = [0 1 3 6 8 10];% change this as desired!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = [0:1/fs:0.01];
pulse = zeros(size(t));
pulse(t>0 & t<=durPulse) = 1;
pulse(t>durPulse & t<=2*durPulse) = -1;
tPulse = [0:1/fs:durPulse];
tPulse = tPulse(1:length(t(t>durPulse & t<2*durPulse)));
[~,endFirst] = min(abs(durPulse-t));
[~,endSecond] = min(abs(2*durPulse -t));

figure
%plot(t,pulse)
hold on

RCVec = [10e-6 10e-5 10e-4 10e-3 10e-2];
modPulse = [];
count = 1;
modPulse = zeros(length(t),length(RCVec));
for RC = RCVec
    modPulse(t>0 & t<=durPulse,count) = pulse(2).*(1-exp(-t(t>0 & t<=durPulse)/RC));
    modPulse(t>durPulse & t<2*durPulse,count) = pulse(endFirst+1)-(pulse(endFirst+1)-modPulse(endFirst-1,count)).*exp(-tPulse/RC);
    modPulse(endSecond:end,count) =  modPulse(endSecond-1,count).*exp(-[0:length([endSecond:length(t)])-1]/(fs*RC));
    plot(t,modPulse(:,count),'linewidth',2)
    count = count + 1;
    
end
title('voltage across capacitor')
leg = cellstr(num2str(RCVec'));
legend(leg)

legend
ylim([-2 2])
%%
figure
tissuePulse = [];
hold on
for count = 1:length(RCVec)
    tissuePulse(:,count) = pulse-modPulse(:,count)';
    plot(t,tissuePulse,'linewidth',2)
end
leg = cellstr(num2str(RCVec'));
legend(leg)
title('Voltage in tissue')
%%
factorDecimate = round(fs/fDesired);
tDecimate = decimate(t,factorDecimate);
modPulseDec = [];
figure
hold on
figure
countPlot = 1;
for index = circshiftIndices
    
    subplot(1,length(circshiftIndices),countPlot)
    hold on
    for count = 1:length(RCVec)
        modPulseDec(:,count) = decimate(circshift(tissuePulse(:,count),index),factorDecimate);
        plot(tDecimate,modPulseDec(:,count),'linewidth',2)
    end
    xlim([0 0.004])
    title({['Shifted ' num2str(index) ' samples'], ['before decimating to ' num2str(fDesired) ' Hz']})
    countPlot = countPlot+1;
end
leg = cellstr(num2str(RCVec'));
legend(leg)