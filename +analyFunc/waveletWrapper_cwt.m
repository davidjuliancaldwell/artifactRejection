function [powerout,f,phaseAngle] = waveletWrapper_cwt(signal,fs,badChans)
% This is a function to run on the response timing data collected by David
% Caldwell and Jeneva Cronin while in the GRIDLab. This uses matlab cwt
% signal:
%   time x channels x trials
% fs:
%   sampling rate in Hz
% badChans:
%   channels to ignore
%
% OUTPUT
% powerout - freq x time x channel x trial

%
% DJC 5-8-2020

%%%%%%%%%%%%%%%%%%%
num_trials = size(signal,3);

if ~exist('badChans','var')
    badChans = [];
end

badChansMask = ones(size(signal,2),1);
badChansMask(badChans) = 0;
badChansMask = logical(badChansMask);
signalTemp = signal(:,badChansMask,:);

totalChans = size(signal,2);

num_channels = size(signalTemp,2);

for i = 1:num_trials
    
    for ii = 1:num_channels
        % compute f and t once
        data_temp = signalTemp(:,ii,i);
        
        if i == 1
            [powerout_temp, f] = cwt(data_temp, fs);
            poweroutTemp(:,:,ii,i) = abs(powerout_temp).^2;
        end
        [powerout_temp] =cwt(data_temp, fs);
        poweroutTemp(:,:,ii,i) = abs(powerout_temp).^2;
    end
    
    
end

%powerout(:,:,badChansMask,:) = poweroutTemp;
%powerout(:,:,~badChansMask,:) = 0;

powerout = poweroutTemp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set phase angle output to be zero

if ~exist('phase_angle','var')
    phaseAngle = [];
end

end