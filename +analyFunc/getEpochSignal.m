%% getEpochSignal.m
%  jdw - 19JUN2011
%
% Changelog:
%   19JUN2011 - originally written
%   11NOV2013 - added cell functionality
%
% This function breaks signal in to epoch based chunks as specified by 
%   starts and ends.  
%
% Parameters:
%   signal - a vector length M containing the signal to break in to epochs.
%   starts - a vector containing the start offsets of each epoch
%   ends - a vector containing the end offsets of each epoch
%
% Return Values:
%   epochSignal - a matrix containing the signal, divided in to epochs.
%     OR
%   epochSignal - a cell array containing the signal, divided in to epochs.
%
% in the case that all of the epochs are not of equal length, this function
% will return a cell array as opposed to a matrix.
%

function epochSignal = getEpochSignal(signal, starts, ends)
%     if (size(signal,2) ~= 1) 
%         warning('2(+)-D matrix passed for signal, ignoring all but first column');
%         while (size(signal, 2) ~= 1)
%             signal = squeeze(signal(:, 1));
%         end
%     end
    if(length(starts) ~= length(ends))
        error('starts and ends must be of same length');
    end
    if (max(ends-starts) ~= min(ends-starts))
        returnMatrix = false;
    else
        returnMatrix = true;
    end

    if (isempty(starts))
        if (returnMatrix == true)
            epochSignal = [];
        else
            epochSignal = {};
        end        
        
        return;
    end
    
    if (returnMatrix == true)
        epochSignal = zeros(ends(1)-starts(1), size(signal,2), length(starts));
    else
        epochSignal = cell(size(signal, 2), length(starts));
    end
    
    for c = 1:length(starts)
        if (returnMatrix == true)
            epochSignal(:,:,c) = signal(starts(c):ends(c)-1,:);
        else
            for chan = 1:size(signal, 2)
                epochSignal{chan, c} = signal(starts(c):ends(c)-1, chan);
            end; clear chan;
        end
    end; clear c;
end