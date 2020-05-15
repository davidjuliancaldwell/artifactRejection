function [] = visualize_wavelet_channel_no_raw_button(powerout,tMorlet,fMorlet,processedSig,tEpoch,chanInt,individual,average,xlims)
% set colormap using cbrewer
%CT = cbrewer('div','RdBu',11);
% flip it so red is increase, blue is down
%CT = flipud(CT);
load('america');
CT = cm;
xlimsVec = xlims;
if individual
    
    for i = 1:size(powerout,4)
        totalFig = figure;
        totalFig.Units = 'inches';
        totalFig.Position = [1 1 4 6];
        subplot(2,1,2);
        s = surf(1e3*tMorlet,fMorlet,powerout(:,:,chanInt,i),'edgecolor','none');
        % Extract X,Y and Z data from surface plot
        x=s.XData;
        y=s.YData;
        z=s.ZData;
        ylimsVec = [5 300];
        %%Create vectors out of surface's XData and YData
        x=x(1,:);
        y=y(1,:);
        interestX = [xlimsVec(1) xlimsVec(end)];
        interestY = [ylimsVec(1) ylimsVec(end)];
        % Divide the lengths by the number of lines needed
        
        % Plot the mesh lines
        % Plotting lines in the X-Z plane
        hold on
        for i = 1:2
            Y1 = interestY(i)*ones(size(x)); % a constant vector
            Z1 = zeros(size(x));
            plot3(x,Y1,Z1,'-k');
        end
        % Plotting lines in the Y-Z plane
        for i = 1:2
            X2 = interestX(i)*ones(size(y)); % a constant vector
            Z1 = zeros(size(X2));
            plot3(X2,y,Z1,'-k');
        end
        
        view(0,90);
        axis tight;
        xlabel('Time (ms)');
        ylabel('Frequency (Hz)');
        title(['Wavelet Decomposition Channel ' num2str(chanInt) ' Trial ' num2str(i)]);
        xlim(xlimsVec);
        ylim(ylimsVec);
        set(gca,'fontsize',8)
        colormap(CT);
        vizFunc.set_colormap_threshold(gcf, [-0.2 0.2], [-2 2], [1 1 1])
        
        h1 = subplot(2,1,1);
        plot(1e3*tEpoch,1e3*processedSig(:,chanInt,i),'color',[204 85 0]/255,'linewidth',1)
        xlabel('Time (ms)');
        ylabel('Voltage (mV)')
        title(['Recovered Signal Channel ' num2str(chanInt) ' Trial ' num2str(i)]);
        ylims = [-(max(abs(1e3*nanmean(squeeze(processedSig(:,chanInt,:)),2)))+0.02) (max(abs(1e3*nanmean(squeeze(processedSig(:,chanInt,:)),2))))+0.02];
        ylim(ylims);
        ylim_h1 = ylims;
        xlim(xlims);
        set(gca,'fontsize',8,'Xlabel',[])
        
    end
    
end
% now average
if average
    
    poweroutAvg = mean(squeeze(powerout(:,:,chanInt,:)),3);
    
    totalFig2 = figure;
    totalFig2.Units = 'inches';
    totalFig2.Position =[1 1 4 6];
    subplot(2,1,2);
    s = surf(1e3*tMorlet,fMorlet,poweroutAvg,'edgecolor','none');
    hold on
    ylimsVec = [5 300];
    % Extract X,Y and Z data from surface plot
    x=s.XData;
    y=s.YData;
    z=s.ZData;
    
    %%Create vectors out of surface's XData and YData
    x=x(1,:);
    y=y(1,:);
    interestX = [xlimsVec(1) xlimsVec(end)];
    interestY = [ylimsVec(1) ylimsVec(end)];
    % Divide the lengths by the number of lines needed
    
    % Plot the mesh lines
    % Plotting lines in the X-Z plane
    hold on
    for i = 1:2
        Y1 = interestY(i)*ones(size(x)); % a constant vector
        Z1 = zeros(size(x));
        plot3(x,Y1,Z1,'-k');
    end
    % Plotting lines in the Y-Z plane
    for i = 1:2
        X2 = interestX(i)*ones(size(y)); % a constant vector
        Z1 = zeros(size(X2));
        plot3(X2,y,Z1,'-k');
    end
    
    view(0,90);
    
    axis tight;
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');
    title(['Time-Frequency Plot Channel ' num2str(chanInt)]);
    xlim(xlimsVec);
    ylim(ylimsVec);
    set(gca,'fontsize',8)
    colormap(CT);
    vizFunc.set_colormap_threshold(gcf, [-0.01 0.01], [-1 1], [1 1 1])
    colorbar();
    
    h3 = subplot(2,1,1);
    plot(1e3*tEpoch,1e3*nanmean(squeeze(processedSig(:,chanInt,:)),2),'color',[204 85 0]/255,'linewidth',1)
    xlabel('Time (ms)');
    ylabel('Voltage (mV)')
    title(['Recovered Signal Channel ' num2str(chanInt)]);
    ylims = [-(max(abs(1e3*nanmean(squeeze(processedSig(:,chanInt,:)),2)))+0.02) (max(abs(1e3*nanmean(squeeze(processedSig(:,chanInt,:)),2))))+0.02];
    ylim(ylims);
    xlim(xlims);
    
    set(gca,'fontsize',8,'Xlabel',[])
    
end

end