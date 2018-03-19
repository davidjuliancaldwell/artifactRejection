function [templateNoEP] = recover_EP(template,varargin)

p = inputParser;

addRequired(p,'template',@isnumeric);
addParameter(p,'plotIt',0,@(x) x==0 || x ==1);

p.parse(template,varargin{:});

template = p.Results.template;
plotIt = p.Results.plotIt;

templateNoEP = template;

for ind = 1:size(template,2)
    templateInd = template(:,ind);
    last = find(abs(zscore(templateInd))>0.005,1,'last');
    last2 = find(abs(diff(templateInd))>20e-6,1,'last')+1;
    ct = max(last,last2);
    
    if length(templateInd)- ct > 10
        x = [ct:length(templateInd)]';
        y = templateInd(x);
        [f2,gof,output] = fit(x,y,'exp2');
        func_fit = @(x) f2.a*exp(f2.b*x) + f2.c*exp(f2.d*x);
        
        % if its a bad fit, set it to be equal to that to try and recover
        % that EP
        
        if gof.adjrsquare<0.5
            templateInd(x) = func_fit(x);
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