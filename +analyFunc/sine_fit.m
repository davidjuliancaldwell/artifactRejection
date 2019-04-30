function [subtractedSig,phase_at_0,f,Rsquare,FITLINE] = sine_fit(X,t,smoothSpan,fRange,fs,plotIt)

% DJC - 1-25-2017 - function takes in a time x channels x trials matrix , and returns
% matrices representing

tRange = 1./fRange;
tRange = tRange.*fs;

if plotIt
 %   figure
end

% make matrices to store values over each iteration
phase_at_0 = zeros(size(X,2),size(X,3),1);
FITLINE = zeros(size(X,1),size(X,2),size(X,3));
subtractedSig = zeros(size(X,1),size(X,2),size(X,3));
Rsquare = zeros(size(X,2),size(X,3),1);
f = zeros(size(X,2),size(X,3),1);

for jj = 1:size(X,2)
    parfor kk = 1:size(X,3)
        sigInd = X(:,jj,kk);
        
        sigIndTemp = sigInd;
        
        maxVal = max(sigInd(zscore(sigInd)<0.5));
        minVal = min(sigInd(zscore(sigInd)>-0.5));
        sigIndTemp(zscore(sigInd)>0.5) = maxVal;
        sigIndTemp(zscore(sigInd)<-0.5) = minVal;
        
        [pha_a,T_a,amp_a,rsquare_a,fitline,offset] = analyFunc.sinfit(1e6*sigIndTemp,smoothSpan,tRange);
        
        f_calculated = 1/(T_a/fs);
        length_sig = length(sigInd);
        x = 1:length_sig;
        a = amp_a.*sin(pha_a+(2*pi*x/T_a))+offset;
        
        phase_delivery = mod((pha_a+(length_sig*pi*2/T_a)),(2*pi));
        
        % save it to matrix
        phase_at_0(jj,kk) = phase_delivery;
        Rsquare(jj,kk) = rsquare_a;
        FITLINE(:,jj,kk) = fitline;
        f(jj,kk) = f_calculated;
        %     if raw_sig
        %         signal_filt = smooth(sig_ind,smooth_span,'moving',0); % smooth data via moving average
        %         signal_filt = signal_filt-median(signal_filt); % DC correction
        %         figure
        %         plot(signal_filt)
        %         hold on
        %         plot(sig_ind)
        %         legend({'smoothed signal','raw signal'})
        %     end
        
        subtractedSig(:,jj,kk) = sigInd' - 1e-6*a;
        
%         if plotIt
%             figure
%             plot(t,a,t,1e6*sigInd,t,1e6*subtractedSig(:,jj,kk))
%             legend({'curve fit','original sig','subtracted'})
%             title('Curve fitting')
%             set(gca,'fontsize',14);
%             %pause(1)
%         end
        fprintf(['Channel ' num2str(jj) ' trial ' num2str(kk) ' complete \n'])
    end
    
end