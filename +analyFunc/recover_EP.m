function [templateNoEP] = recover_EP(template,fs,varargin)

p = inputParser;

addRequired(p,'template',@isnumeric);
addRequired(p,'fs',@isnumeric);
addParameter(p,'plotIt',0,@(x) x==0 || x ==1);


p.parse(template,fs,varargin{:});

template = p.Results.template;
fs = p.Results.fs;
plotIt = p.Results.plotIt;

templateNoEP = template;
pctl = @(v,p) interp1(linspace(0.5/length(v), 1-0.5/length(v), length(v))', sort(v), p*0.01, 'spline');
timeSampsExtend = 1*fs/1000;% time_ms

for ind = 1:size(template,2)
    ct = [];
    templateInd = template(:,ind);
    
    absZSig = abs(zscore(templateInd));
    absZDiffSig = abs(zscore(diff(templateInd)));
    threshSig = pctl(absZSig,97.5);
    threshDiff = pctl(absZDiffSig,97.5);
    
    last = find(absZSig>threshSig,1,'last'); % started with 0.2
    last2 = find(absZDiffSig>threshDiff,1,'last')+1; % started with 5e-3
    
    
    %     last = find(abs(zscore(templateInd))>0.005,1,'last');
    %     last2 = find(abs(diff(templateInd))>20e-6,1,'last')+1;
    ct = max(last,last2);
    
    if (length(templateInd)- ct > timeSampsExtend) & ~isempty(ct)
        x = [ct:length(templateInd)]';
        y = templateInd(x);
        try
            [f2,gof,output] = fit(x,y,'exp2');
            func_fit = @(x) f2.a*exp(f2.b*x) + f2.c*exp(f2.d*x);
            
            % if its a bad fit, set it to be equal to that to try and recover
            % that EP
            
            if gof.adjrsquare<0.5
                templateInd(x) = func_fit(x);
            end
            
        catch
            templateInd = templateInd;
            
        end
        
        if plotIt
            figure
            plot(templateInd)
            hold on
            plot(x,y,'linewidth',2)
            plot(x,func_fit(x));
            vline(last2)
            vline(last,'g')
            
        end
    end
    templateNoEP(:,ind) = templateInd;
    
    
end

end