function plot_error(x, y,type,color,edge,add,transparency)

% plot standard error around mean
% assume x is time, y is time x trials

if nargin<7;transparency=.3;end %default is to have a transparency of .3
if nargin<6;add=1;end     %default is to add to current plot
if nargin<5;edge='k';end  %dfault edge color is black
if nargin<4;color='b';end %default color is blue
if nargin<3;type='SE';end %default is standard error, otherwise standard deviation

washeld = ishold(gca);

sigma = std(y,0,2);
n = size(y,2);
if strcmp(type,'SE')% standard error
    upper = mean(y,2) + sigma/sqrt(n);
    lower = mean(y,2) - sigma/sqrt(n);
elseif strcmp(type,'SD') % standard deviation
    upper = mean(y,2) + sigma;
    lower = mean(y,2) - sigma;
elseif strcmp(type,'CI')
    upper = mean(y,2) + 1.96*sigma/sqrt(n); % 95% confidence interval
    lower = mean(y,2) - 1.96*sigma/sqrt(n);
end

if (washeld == false)
    hold on;
end

%plot(x, squeeze(nanmean(y,2)), color);
plot(x, squeeze(nanmean(y,2)), 'color',color) % djc change 2-5-2018

[fillhandle,msg]=vizFunc.jbfill(x,upper',lower',color',edge,add,transparency);

if (washeld == false)
    hold off;
end

axis tight;

end
