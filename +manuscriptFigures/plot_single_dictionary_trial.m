%% script to plot dictionary extracted vs. all the trials
%
% David.J.Caldwell 9.21.2018

chanInt = 47;
figure
hold on
timeVec = 1e3*[0:size(templateTrial{chanInt}{:,i},1)-1]/fsData;

plot(timeVec,templateDictCell{chanInt},'linewidth',8,'color','k')
for i = 1:size(templateTrial{chanInt},2)
    plot(timeVec,templateTrial{chanInt}{:,i})
end
xlim([0 2])
title(['Channel ' num2str(chanInt)])
title({'Extracted templates',' and individual trial artifacts'})
xlabel('time (ms)')
ylabel('Voltage (V)')
set(gca,'fontsize',20)

%%

chanInt = 32;
figure
for i = 1:size(templateTrial{chanInt},2)
    timeVec = 1e3*[0:size(templateTrial{chanInt}{:,i},1)-1]/fs;
    
    vizFunc.plot_error(timeVec,templateTrial{chanInt}{:,i},'CI',rand(1,3));
end
title(['Channel ' num2str(chanInt)])
xlabel('time (ms')
ylabel('Voltage (V)')
