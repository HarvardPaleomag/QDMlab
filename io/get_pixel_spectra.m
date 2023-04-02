function [lpos, rpos, lneg, rneg, freqList] = get_pixel_spectra(rawDataPos, rawDataNeg, idx)
% [lpos, rpos, lneg, rneg, freqList] = get_pixel_spectra(rawDataPos, rawDataNeg, idx)
% Extract left and right spectra for positive and negative data sets at the specified index.
% 
% Parameters:
% ----------
% - rawDataPos: struct containing the raw data for the positive set, with fields:
%   * imgStack1: 3D matrix of image data for left camera
%   * imgStack2: 3D matrix of image data for right camera
%   * numFreqs: number of frequency points in the image data
%   * imgNumCols: number of columns in the image data
%   * imgNumRows: number of rows in the image data
%   * freqList: 1D vector of frequency values for the image data
% - rawDataNeg: struct containing the raw data for the negative set, with fields:
%   * imgStack1: 3D matrix of image data for left camera
%   * imgStack2: 3D matrix of image data for right camera
%   * numFreqs: number of frequency points in the image data
%   * imgNumCols: number of columns in the image data
%   * imgNumRows: number of rows in the image data
%   * freqList: 1D vector of frequency values for the image data
% - idx: 1x2 vector specifying the index of the pixel to extract spectra from
%
% Returns:
% --------
% - lpos: column vector representing the left spectrum for the positive set at the specified index
% - rpos: column vector representing the right spectrum for the positive set at the specified index
% - lneg: column vector representing the left spectrum for the negative set at the specified index
% - rneg: column vector representing the right spectrum for the negative set at the specified index
% - freqList: 1D vector representing the frequency values for the image data

    % Reshape the image data from both sets into a 3D matrix with dimensions (numFreqs, imgNumCols, imgNumRows)
    dataPosLeft = reshapeImg(rawDataPos.imgStack1, rawDataPos.numFreqs, ...
        rawDataPos.imgNumCols, rawDataPos.imgNumRows);
    dataPosRight = reshapeImg(rawDataPos.imgStack2, rawDataPos.numFreqs, ...
        rawDataPos.imgNumCols, rawDataPos.imgNumRows);
    dataNegLeft = reshapeImg(rawDataNeg.imgStack1, rawDataNeg.numFreqs, ...
        rawDataNeg.imgNumCols, rawDataNeg.imgNumRows);
    dataNegRight = reshapeImg(rawDataNeg.imgStack2, rawDataNeg.numFreqs, ...
        rawDataNeg.imgNumCols, rawDataNeg.imgNumRows);
    
    % Extract the frequency list from the positive set
    freqList = reshape(rawDataPos.freqList, [rawDataPos.numFreqs, 2]);

    % Extract the left and right spectra for both the positive and negative sets at the specified index
    lneg = dataNegLeft(:, idx(1), idx(2));
    rneg = dataNegRight(:, idx(1), idx(2));
    lpos = dataPosLeft(:, idx(1), idx(2));
    rpos = dataPosRight(:, idx(1), idx(2));
end
