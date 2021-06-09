function [corrected, debug] = correct_global_vec(data, globalFraction, kwargs)
%[corrected, debug] = correct_global(data, globalFraction; 'mean')
% global spectra subtraction
% Find global averaged spectra

arguments
    data
    globalFraction
    kwargs.mean = 'none'
end

if strcmp(kwargs.mean, 'none')
    specData = squeeze(mean(data,2));%global spectrum 51x1 array
else
    specData = kwargs.mean;
end

% BL = base line
leftBL = mean(specData(1:5,1));
rightBL = mean(specData(end-5:end,1));

globalmeanBL = mean([leftBL rightBL]);
baselinerange = globalmeanBL-min(specData);%difference between baseline and the minimum value in global spectrum

specZBL = specData-globalmeanBL; %ZBL=zero baseline

corrected = zeros(size(data));

debug = zeros([size(data,2),1]);

leftBL = mean(data(1:5,:),1);
rightBL = mean(data(size(specData,1)-5:end,:),1);
meanBL = mean([leftBL; rightBL]);

pixelrange = meanBL-min(data); %difference between baseline and the minimum value in pixel spectrum

baselineRatio = (baselinerange ./ pixelrange);
corrected = baselineRatio .* (data - meanBL) - globalFraction*specZBL + globalmeanBL;
