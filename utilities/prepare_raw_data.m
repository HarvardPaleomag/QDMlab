function [binDataNorm, freq] = prepare_raw_data(expData, binSize, nRes, kwargs)
%[binDataNorm, freq] = prepare_raw_data(expData, binSize, nRes; 'gpuData')
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

arguments
    expData
    binSize
    nRes
    kwargs.gpuData = false
end

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
    if nRes == 1
        a = 1; b = 2;
    else
        a = 3; b = 4;
    end
    dataStacka = expData.(sprintf('imgStack%i',a)); 
    dataStackb = expData.(sprintf('imgStack%i',b));
    dataStack = [dataStacka; dataStackb];
end

% reshape and transpose each image
% NOTE: I guess the data is written pixel by pixel, starting left top
% then row by row -> freq x col * rows (i.e. [51, 2304000])
data = reshape(dataStack, [], expData.imgNumCols, expData.imgNumRows); % reshape into freq x col x rows
data = permute(data,[3 2 1]); % permute the axis to rows x cols x freq

% binning
msg = sprintf('<>   %i: binning data >> binSize = %i', nRes, binSize);
logMsg('debug',msg,1,0);

binData = imresize(data, 1/binSize, 'method', 'box');

% Correct for severely non-unity baseline by dividing pixelwise by
% average of all frequency points

binDataNorm = zeros(size(binData));
NormalizationFactor = mean(binData,3);    % compute average

for y = 1:length(freq)
    binDataNorm(:,:,y) = binData(:,:,y) ./ NormalizationFactor;
end

%% return gpudata if kwargs.gpuData == true
if isequal(kwargs.gpuData, true)
    [sizeY, sizeX, sweepLength] = size(binDataNorm, [1,2,3]);
    imgPts = sizeX*sizeY;
    gpudata = reshape(binDataNorm, [imgPts, sweepLength]); % make it into 2d matrix
    gpudata = transpose(gpudata); %transpose to make it 51 x pixels
%     gpudata = single(gpudata);
    binDataNorm = gpudata;
end

