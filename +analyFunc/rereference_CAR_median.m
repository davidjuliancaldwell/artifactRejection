function output=rereference_CAR_median(data,mode,badChannels,permuteOrder,channelReference,channelsToUse)
% Function to either common average, or median rereference,
% while excluding bad channels.
%
% data is either time x channels, or time x channels x trials
% bad_channels is a list of bad channels as output
% leaves output channels that are "bad" untouched
% permute is order if need to permute, like [1 3 2]
% channelsToUse are the ones to use if
%
% David. J.Caldwell 6.12.2017

if (~exist('badChannels','var') || isempty(badChannels))
    badChannels = []; % default no bad channels, this overrides good channels
end

if (~exist('permuteOrder','var') || isempty(permuteOrder))
    permuteOrder = [1:ndims(data)];
end

if (~exist('channelReference','var')|| isempty(channelReference))
    channelReference = [1]; % default to rereference against single channel if need be
end

if (~exist('channelsToUse','var') || isempty(channelsToUse))
    channelsToUse = logical(ones(size(data,2),1));
end

data = permute(data,permuteOrder);

output = data; % default output of data with permute order

channelMask = logical(ones(size(data,2),1));
channelMask(badChannels) = 0;

switch(mode)
    case 'mean'
        avg = mean(data(:,channelMask,:),2);
        avg = repmat(avg, 1, size(data(:,channelMask,:),2));
        output(:,channelMask,:) = data(:,channelMask,:) - avg;
        
        % shift data if needed
        output = permute(output,permuteOrder);
    case 'median'
        med = median(data(:,channelMask,:),2);
        med = repmat(med, 1, size(data(:,channelMask,:),2));
        output(:,channelMask,:) = data(:,channelMask,:) - med;
        
        % shift data if needed
        output = permute(output,permuteOrder);
        
    case 'singleChan'
        chan = data(:,channelReference,:);
        repmatChan = repmat(chan,1,size(data,2));
        output = data - repmatChan;
        
        output = permute(output,permuteOrder);
        
        
    case 'bipolarPair' % do 1 vs 2, 3 vs 4, etc
        for i = [1:2:size(data,2)]
            chanOdd = data(:,i,:);
            chanEven = data(:,i+1,:);     
            newChan = chanEven-chanOdd;
            output(:,i,:) = zeros(size(newChan));
            output(:,i+1,:) = newChan;   
        end
        
        output = permute(output,permuteOrder);
        
        
    case 'bipolar' % do 1 vs 2, 3 vs 4, etc
        for i = [1:7] % for R side , do 8:end
            chanOdd = data(:,i,:);
            chanEven = data(:,i+1,:);        
            newChan = chanEven-chanOdd;
            output(:,i,:) = newChan;
        end
        output(:,8,:) = zeros(size(newChan));
        output = permute(output,permuteOrder);
    case 'selectedChannels'
        avg = mean(data(:,channelsToUse,:),2);
        avg = repmat(avg, 1, size(data(:,channelMask,:),2));
        output(:,channelMask,:) = data(:,channelMask,:) - avg;
        output = permute(output,permuteOrder);
        
    case 'none'
        output = data;
        
    case ''
        output = data;
        
end