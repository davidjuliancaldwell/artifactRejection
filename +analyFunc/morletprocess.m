function [ powerout, f, t, phaseangle ] = morletprocess( inputs, fs, time_res, use_lmp )
%MORLETPROCESS Processing ECoG spectrum using analytic/nonanalytic Morlet
%wavelets in Matlab cwtft(). Outputs spectral power.
%   [POWEROUT, F, T, PHASEANGLE] = MORLETPROCESS(INPUTS, FS, TIME_RES, USE_LMP)
%
%   This function wraps cwtft to get an efficiently spaced Morlet/Morlex
%   power spectrum estimation with certain time resolution bins.
%
%   F (frequency) and T (time) bin labels are provided in the output.
%
%   Using this function to compute phase is not necessarily recommended
%   as Hilbert phases are generally considered more stable for broadband
%   behavior, but can be used if needed. cwtft 'morlex' (non-analytic) is
%   used if phase is not required, otherwise 'morl' is used.
%
%   POWEROUT = MORLETPROCESS(INPUTS, FS, TIME_RES) takes INPUTS as
%   [time x channels] with time along dim 1. It iterates through all
%   channels to compute POWEROUT, with sampling frequency FS, and bins by
%   TIME_RES in seconds.
%
%   Example: POWEROUT = MORLETPROCESS(INPUTS, 1200, 0.050) for 1200 Hz
%   sampling frequency and 50ms bins.
%
%   USE_LMP is a bool (true, false). Setting USE_LMP to true will add an
%   additional "0 Hz" frequency track to all power outputs, to estimate
%   "local motor potential".
%

    if(~exist('use_lmp', 'var'))
        use_lmp = false; % do not add LMP "0 Hz" track by default
    end
    
    if(nargout == 4)
        use_wavelet = 'morl'; % only if phase angle is required
    else
        use_wavelet = 'morlex'; % non-analytic, runs faster
    end
    
    % initialize dt and bins
    dt = 1/fs;
    binsize = round(fs*time_res);
    truncateby = mod(length(inputs), binsize); % make round bins
    
    % establish scales, frequencies, and time
    numvoices = 4;
    a0 = 2^(1/numvoices);
    f0 = 3/pi;
    minfreq = 1;
    maxfreq = 300;

    minscale = f0/(maxfreq/fs);
    maxscale = f0/(minfreq/fs);
    minscale = floor(numvoices*log2(minscale));
    maxscale = ceil(numvoices*log2(maxscale));
    scales = a0.^(minscale:maxscale)./fs;

    f = f0./scales;
    t = (binsize*dt/2):(binsize*dt):((length(inputs)-truncateby)*dt);

    % initialize power, phase, and lmp tracks, if any
    powerout = zeros(length(f), length(t), size(inputs, 2));
    if(nargout == 4)
        phaseangle = zeros(size(powerout));
    end
    if(use_lmp)
        siglmp = zeros(1, length(t), size(inputs, 2));
    end
    
    for i = 1:size(inputs, 2)
        cwty = cwtft({inputs(:, i), dt},'wavelet',use_wavelet,'scales',scales,'padmode','symw');
        powerout(:, :, i) = truncbindata(abs(cwty.cfs).^2, truncateby, binsize);
        phaseangle(:, :, i) = truncbindata(angle(cwty.cfs), truncateby, binsize);

        if(use_lmp)
            lmpprebin = computeLMP(inputs(:, i), fs);
            siglmp(:, :, i) = truncbindata(lmpprebin, truncateby, binsize);
        end
    end

    if(use_lmp) % stick LMP on to the end of power, as a "0 Hz" track
        powerout = cat(1,powerout, siglmp);
        f(end+1) = 0;
    end
end

function [ postbin ] = truncbindata( prebin, truncateby, binsize, forcedim )
% This function bins input along time, and assumes time is the longer dim.
% If this is not true, set forcedim = 1 or 2
%
% When output, time is along dim 2
    if(~exist('forcedim', 'var'))
        forcedim = [];
    end
    if(isempty(forcedim))
        if(size(prebin, 2) > size(prebin, 1))
            prebin = prebin';
        end
    elseif(forcedim == 2)
        prebin = prebin';
    elseif(forcedim ~= 1)
        error('dimensional error on binning');
    else
        error('dimensional parameter error on binning');
    end

    binningtemp = prebin(1:end-truncateby, :);
    binningtemp = reshape(binningtemp, binsize, size(binningtemp, 1)/binsize, size(binningtemp, 2));
    binningtemp = squeeze(mean(binningtemp, 1));
    
    postbin = binningtemp';
end

function [ out_smooth ] = computeLMP( in_smooth, fs )
% Uses specifically designed lowpass to obtain LMP. This obtains the
% frequency response of moving average windows of ~133ms (Schalk 2007)
% using a low-order zero-phase Butterworth filter to achieve minimum phase
% lags and smooth frequency response characteristics.

    v = version;
    if(str2double(v(1:3)) >= 8.3)
        
    out_smooth = zeros(size(in_smooth));
    
    % MATLAB 2014a and after
        d = designfilt('lowpassiir', 'PassbandFrequency', 0.2, ...
                                    'StopbandFrequency', 1/2, ...
                                    'PassbandRipple', 1, ...
                                    'StopbandAttenuation', 3, ...
                                    'SampleRate', fs, ...
                                    'DesignMethod', 'butter');

        for i = 1:size(in_smooth, 2)
            out_smooth(:, i) = filtfilt(d, in_smooth(:, i));
        end
        
    else

    % MATLAB 2013b and before
        
        h  = fdesign.lowpass(0.2, 1/2, 1, 3, fs);
        Hd = design(h, 'butter', 'MatchExactly', 'stopband');
        [b, a] = sos2tf(Hd.sosMatrix, Hd.ScaleValues);

        for i = 1:size(in_smooth, 2)
            out_smooth(:, i) = filtfilt(b, a, in_smooth(:, i));
        end
        
    end
end
