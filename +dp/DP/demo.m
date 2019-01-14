
% demos 
 
%% Collapse Gibbs sampling for Dirichelt process gaussian mixture model
%close all; clear;
d = 2;
k = 3;
n = 500;
[XX,label] = mixGaussRnd(d,k,n);
plotClass(XX,label);

[y,model] = mixGaussGb(XX);
figure
plotClass(XX,y);

% 
