function [binDataNorm, freq] = prepare_raw_data(expData, binSize, nRes)
% [binDataNorm, freq] = prepare_raw_data(expData, binSize, nRes)
% prepares the raw data for GPU fitting
% 1. reshapes the data into from (x*y) -> (y,x) array
% 2. Bins data: (imresize)
% 3. Normalizes
%
% Parameters
% ----------
%     required
%     ========
%     expData: struct
%         Data of load(run0000n.mat)
%     binSize: int
%         binning size (can be 1)
%     nRes: int
%         number of resonance. Low frequencies = 1, High frequencies = 2

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
binData = imresize(data, 1/binSize, 'method', 'box');

% Correct for severely non-unity baseline by dividing pixelwise by
% average of all frequency points

binDataNorm = zeros(size(binData));
NormalizationFactor = mean(binData,3);    % compute average

for y = 1:length(freq)
    binDataNorm(:,:,y) = binData(:,:,y) ./ NormalizationFactor;
end
