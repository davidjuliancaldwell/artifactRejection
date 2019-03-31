function [] = visualize_raw_vs_processed(processedSig,tEpoch,dataInt,chanInt,individual,average)

if individual
    
    for i = 1:size(processedSig,3)
 
        h1 = subplot(2,1,2);
        plot(1e3*tEpoch,1e3*processedSig(:,chanInt,i),'color',[204 85 0]/255,'linewidth',0.5)
        xlabel('Time (ms)');
        ylabel('Voltage (mV)')
        title(['Recovered Signal Channel ' num2str(chanInt) ' Trial ' num2str(i)]);
        ylims = [-(max(abs(1e3*processedSig(:,chanInt,i))) + 10) (max(abs(1e3*processedSig(:,chanInt,i))) + 10)];
        ylim(ylims);
        ylim_h1 = ylims;
        xlim([-200 1000]);
        set(gca,'fontsize',14,'Xlabel',[])
        
        h2 = subplot(2,1,1);
        plot(1e3*tEpoch,1e3*dataInt(:,chanInt,i),'color','k','linewidth',3)
        xlabel('Time (ms)');
        ylabel('Voltage (mV)')
        title(['Raw Signal Channel ' num2str(chanInt) ' Trial ' num2str(i)]);
        ylim(ylim_h1);
        xlim([-200 1000]);
        set(gca,'fontsize',14,'Xlabel',[])
        
        linkaxes([h1,h2],'xy');
        
    end
    
end
% now average
if average
       
    totalFig2 = figure;
    totalFig2.Units = 'inches';
    totalFig2.Position = [1 1 4 4];
        
    h3 = subplot(2,1,2);
    plot(1e3*tEpoch,1e3*nanmean(squeeze(processedSig(:,chanInt,:)),2),'color',[204 85 0]/255,'linewidth',0.5)
    xlabel('Time (ms)');
    ylabel('Voltage (mV)')
    title(['Recovered Signal Channel ' num2str(chanInt)]);
    ylims = [-(max(abs(1e3*nanmean(squeeze(processedSig(:,chanInt,:)),2))) + 10) (max(abs(1e3*nanmean(squeeze(processedSig(:,chanInt,:)),2))) + 10)];
    ylim(ylims);
    ylim_h1 = ylims;
    xlim([-200 1000]);
    
    set(gca,'fontsize',12,'Xlabel',[])
    
    h4 = subplot(2,1,1);
    plot(1e3*tEpoch,1e3*nanmean(squeeze(dataInt(:,chanInt,:)),2),'color','k','linewidth',0.5)
    xlabel('Time (ms)');
    ylabel('Voltage (mV)')
    title(['Raw Signal Channel ' num2str(chanInt)]);
    ylim(ylim_h1);
    xlim([-200 1000]);
    set(gca,'fontsize',12,'Xlabel',[])
    
    linkaxes([h3,h4],'xy');
end

end