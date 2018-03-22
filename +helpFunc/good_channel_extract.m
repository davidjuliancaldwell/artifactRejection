function [goods,goodVec] = good_channel_extract(varargin)
% Usage:  [goods,goodVec] = good_channel_extract(varargin)
%
% This function will perform an interpolation scheme for artifacts on a
% trial by trial, channel by channel basis, implementing either a linear
% interpolation scheme, or a pchip interpolation scheme
% 
% Arguments:
%   Required:
%        numChans - Number of channels recorded from, including stimulation
%        and bad channels
%
%   Optional:
%   bads - samples x channels x trials 
%
%         stimChans - channels used for stimulation which should be
% Returns:
%   goods - a logical matrix with 0 being the "bad" channels, and 1 being
%   good channels (e.g. [1 0 0 1]
%   goodVec - vector of the channel indices that are good (e.g. [4 10 15
%   20])
%
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

p = inputParser;

addParameter(p,'numChans',64,@isnumeric);
addParameter(p,'bads',[],@isnumeric);
addParameter(p,'stimChans',[],@isnumeric);
p.parse(varargin{:});

bads = p.Results.bads;
stimChans = p.Results.stimChans;
numChans = p.Results.numChans;

badTotal = [stimChans; bads];

% make logical good channels matrix to index
goods = zeros(numChans,1);
channelsOfInt = 1:numChans;
goods(channelsOfInt) = 1;

% set the goods matrix to be zero where bad channels are
goods(badTotal) = 0;

% make it logical
goods = logical(goods);
[~,goodVec] = find(goods'>0);

end

