%% ecogFilter.m
%  jdw - 28APR2011
%
% Changelog:
%   28APR2011 - originally written
%
% This function filters an input signal with multiple notch filters and a
% high pass and low pass filter if desired.  It can do each of these
% independent of the others.
%
% Parameters:
%   signals - signals to be filtered.  If signals is an MxN array, the
%     vectors of length M, indexed by the N dimension will be treated as
%     independent signals and filtered as such.
%   lnReject - a boolean value that, if true indicates that the signals
%     should be filtered with notch filters centered around lnFreqs
%   lnFreqs - the frequencies that will be used for notch filters if
%     lnReject is true
%   hp - a boolean value that, if true indicates that the signals should be
%     filtered with a high pass filter
%   hpFreq - the highpass cutoff frequency
%   lp - a boolean value that, if true indicates that the signals should be
%     filtered with a low pass filter
%   lpFreq - the lowpass cutoff frequency
%   fSamp - the sampling rate of signals
%   filterOrder (optional) - the filter order of the butterworth filter to
%     be used.  The default value is 4th order.
%
% Return Values:
%   filteredSignals - the filtered signals
%
function filteredSignals = ecogFilter(signals, lnReject, lnFreqs, hp, hpFreq, lp, lpFreq, fSamp, filterOrder, causality)
    if(~exist('filterOrder','var') || isempty(filterOrder))
        filterOrder = 4;        
    end
    
    if(~exist('causality', 'var'))
        causality = 'acausal';
    end
    
    
    
    if (lnReject == true)
        start_frqs = zeros (size(lnFreqs));
        stop_frqs  = zeros (size(lnFreqs));

        for j=1:size(lnFreqs,2)
          start_frqs(j) = (lnFreqs(j) - 3) / (fSamp/2);     
          stop_frqs(j)  = (lnFreqs(j) + 5) / (fSamp/2);     
          
          [b, a] = butter (filterOrder, [start_frqs(j) stop_frqs(j)], 'stop');

          signals = subFilter(b, a, signals, causality);
        end
    end

    if (hp == true)
        [b, a] = butter (filterOrder, hpFreq / (fSamp/2), 'high');
        signals = subFilter(b, a, signals, causality);
    end
    
    if (lp == true)
        [b, a] = butter (filterOrder, lpFreq / (fSamp/2), 'low');
        signals = subFilter(b, a, signals, causality);
    end
    
    filteredSignals = signals;
end

function rSignals = subFilter(b, a, signals, causality)
    rSignals = zeros(size(signals));
    
    switch(causality)
        case 'acausal'
            rSignals = filtfilt(b, a, signals);
        case 'causal'
           rSignals = filter(b,a,signals); 
    end
end