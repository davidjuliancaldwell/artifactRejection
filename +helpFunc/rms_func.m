function [rms] = rms_func(data)
% takes in data that is time x channels and returns the RMS along the
% dimension of a channel. This may be a useful metric of signal quality 

rms = sqrt(sum(data.^2,1)/size(data,1));

end