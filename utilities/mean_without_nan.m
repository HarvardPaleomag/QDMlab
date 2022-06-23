function meanD = mean_without_nan(data)
%[meanD] = mean_without_nan(data)
% takes a QDM image stack and returns the mean for only those pixels, tat
% do not contain NaN in the data

% find pixels with nan
nanPix = isnan(mean(data,1));

nPix = numel(nonzeros(nanPix));

if nPix >0
    msg = sprintf('Calculating mean without << %i NaN >> pixels.', nPix);
    logMsg('info',msg,1,0);
end

meanD = data(:,~nanPix);
% calculate mean for remaining pixels
meanD = mean(meanD, 2);