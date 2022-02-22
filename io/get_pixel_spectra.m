function [lpos, rpos, lneg, rneg, freqList] = get_pixel_spectra(rawDataPos, rawDataNeg, idx)
%[lpos, rpos, lneg, rneg, freqList] = get_pixel_spectra(rawDataPos, rawDataNeg, idx)
    dataPosLeft = reshapeImg(rawDataPos.imgStack1, rawDataPos.numFreqs, ...
        rawDataPos.imgNumCols, rawDataPos.imgNumRows);
    dataPosRight = reshapeImg(rawDataPos.imgStack2, rawDataPos.numFreqs, ...
        rawDataPos.imgNumCols, rawDataPos.imgNumRows);
    dataNegLeft = reshapeImg(rawDataNeg.imgStack1, rawDataNeg.numFreqs, ...
        rawDataNeg.imgNumCols, rawDataNeg.imgNumRows);
    dataNegRight = reshapeImg(rawDataNeg.imgStack2, rawDataNeg.numFreqs, ...
        rawDataNeg.imgNumCols, rawDataNeg.imgNumRows);
    freqList = reshape(rawDataPos.freqList, [rawDataPos.numFreqs, 2]);

    lneg = dataNegLeft(:, idx(1), idx(2));
    rneg = dataNegRight(:, idx(1), idx(2));
    lpos = dataPosLeft(:, idx(1), idx(2));
    rpos = dataPosRight(:, idx(1), idx(2));
    
