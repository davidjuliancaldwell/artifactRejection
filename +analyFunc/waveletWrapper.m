function [powerout,f,t,phaseAngle] = waveletWrapper(signal,fs,timeRes,badChans)
% This is a function to run on the response timing data collected by David
% Caldwell and Jeneva Cronin while in the GRIDLab. This uses a
% morletprocess script as implemented by James Wu
% signal:
%   time x channels x trials
% fs:
%   sampling rate in Hz
% time_res:
%   time resolution for morlet process
% badChans:
%   channels to ignore
%
% OUTPUT 
% powerout - freq x time x channel x trial

%
% DJC 3-30-2017

% Morlet Process
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

for i = 1:num_trials
    
    data_temp = signalTemp(:,:,i);
    % compute f and t once
    if i == 1
        [powerout_temp, f, t] = analyFunc.morletprocess(data_temp, fs, timeRes);
        poweroutTemp(:,:,:,i) = powerout_temp;
    end
    [powerout_temp] = analyFunc.morletprocess( data_temp, fs, timeRes);
    poweroutTemp(:,:,:,i) = powerout_temp;
    
   
end

powerout(:,:,badChansMask,:) = poweroutTemp;
powerout(:,:,~badChansMask,:) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set phase angle output to be zero 

if ~exist('phase_angle','var')
   phaseAngle = []; 
end

end