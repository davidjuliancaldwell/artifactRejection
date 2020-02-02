function [buttonLocsSampsCell,buttonLocsCell,data,tEpoch,uniqueCond] = combine_pooled_data(sid,epochedCortEco_cell,tEpoch)

block = [1,2];
buttonLocsSampsCellInd = {};
buttonlocsSamps_cell = {};
data = {};
tEpochGood = tEpoch;

for ii = block
    % load([sid,'_compareResponse_block_',num2str(i),'_changePts_tactorSub .mat'])
    load([sid,'_compareResponse_block_',num2str(ii),'_changePts_noDelay.mat'])
    
    buttonLocsSampsCellInd{ii} = buttonLocsSamps; % samples
    buttonLocsCellInd{ii} = buttonLocs; % seconds
end

numConditions = length(buttonLocs);

for jj = 1:numConditions
    buttonLocsSampsCell{jj} = [[buttonLocsSampsCellInd{1}{jj}] [buttonLocsSampsCellInd{2}{jj}]];
    buttonLocsCell{jj} =  [[buttonLocsCellInd{1}{jj}] [buttonLocsCellInd{2}{jj}]];
    data{jj} = cat(3,[epochedCortEco_cell{1}{jj}], [epochedCortEco_cell{2}{jj}]);
    
end

tEpoch = tEpochGood;

end