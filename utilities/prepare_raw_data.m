function [binDataNorm, freq] = prepare_raw_data(expData, binSize, nRes)
% prepares the raw data. Does normalization, binning and reshaping

dataStack = expData.(sprintf('imgStack%i',nRes));

if nRes == 1
    freq = expData.freqList(1:expData.numFreqs) / 1E9;   %everything is in GHz

else
    freq = expData.freqList(1+expData.numFreqs : 2*expData.numFreqs) / 1E9;
end

%% data preparation
% X/Y of unbinned data
% Note: X = COL; Y = ROW -> (y,x) for (row, col) matlab convention
spanXTrans = 1:expData.imgNumCols;
spanYTrans = 1:expData.imgNumRows;

% check for 101 frequencies. File includes imgStack3
if isfield(expData, 'imgStack3')      
    % combine 1&2 or 3&4
    dataStacka = expData.(sprintf('imgStack%i',nRes)); 
    dataStackb = expData.(sprintf('imgStack%i',nRes+1));
    dataStack = [dataStacka; dataStackb];
end

data = zeros(expData.imgNumRows, expData.imgNumCols, expData.numFreqs);

% crop
data = data(spanYTrans,spanXTrans,:);

% reshape and transpose each image
for y = 1:expData.numFreqs
    data(:,:,y) = transpose(reshape(dataStack(y, :), [expData.imgNumCols, expData.imgNumRows]));
end

% binning
fprintf('<>   %i: binning data >> binSize = %i\n', nRes, binSize);

sizeXY = size(BinImage(data(:,:,1),binSize));
binData = zeros(sizeXY(1),sizeXY(2),length(freq));

for y = 1:length(freq)
    binData(:,:,y) = BinImage(data(:,:,y),binSize);
end

sizeX = size(binData,2); % binned image x-dimensions
sizeY = size(binData,1); % binned image y-dimensions

% Correct for severely non-unity baseline by dividing pixelwise by
% average of all frequency points

binDataNorm = zeros(size(binData));
NormalizationFactor = mean(binData,3);    % compute average

for y = 1:length(freq)
    binDataNorm(:,:,y) = binData(:,:,y) ./ NormalizationFactor;
end