function gpudata = gpu_data_reshape(binDataNorm)
%[gpudata] = gpu_data_reshape(binDataNorm)
% function takes the binDataNorm style (e.g. 1200x1920) data array and
% transforms it into the gpu formatted array

[sizeY,sizeX, sweepLength] = size(binDataNorm, [1,2,3]);
imgPts = sizeX*sizeY;
gpudata = reshape(binDataNorm, [imgPts, sweepLength]); % make it into 2d matrix
gpudata = transpose(gpudata); %transpose to make it 51 x pixels
gpudata = single(gpudata);