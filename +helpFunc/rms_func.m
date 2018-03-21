function [rms] = rms_func(data)
% takes in data that is samples x channels and returns the RMS along the
% dimension of a channel. This may be a useful metric of signal quality 
%
% Arguments:
%   Required:
%   data - samples x channels  
%
%
% Returns:
%   rms - rms values for each of the channels
%
%
% Copyright (c) 2018 Updated by David Caldwell
% University of Washington
% djcald at uw . edu 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rms = sqrt(sum(data.^2,1)/size(data,1));

end