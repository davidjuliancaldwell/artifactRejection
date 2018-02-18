function [goods,goodVec] = goodChannel_extract(varargin)
%GOODCHANNEL_EXTRACT This function extracts 

stimChans = [];
bads = [];
numChans = 64;

for i=1:2:(length(varargin)-1)
    switch lower(varargin{i})
        case 'bads'
            bads = varargin{i+1};
        case 'stimchans'
            stimChans = varargin{i+1};
        case 'numchans'
            numChans = varargin{i+1};
    end
end

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

