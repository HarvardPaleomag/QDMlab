function [corrected, debug] = correct_global(data, globalFraction, kwargs)
%[corrected, debug] = correct_global(data, globalFraction; 'mean')
% global spectra subtraction
% Find global averaged spectra

arguments
    data
    globalFraction
    kwargs.mean = 'none'
end

if strcmp(kwargs.mean, 'none')
    specData = squeeze(mean(data,[1 2]));%global spectrum 51x1 array
else
    specData = kwargs.mean;
end

% BL = base line
leftBL = mean(specData(1:5,1));
rightBL = mean(specData(size(specData,1)-5:end,1));
globalmeanBL = mean([leftBL rightBL]);
baselinerange = globalmeanBL-min(specData);%difference between baseline and the minimum value in global spectrum
specZBL = specData-globalmeanBL; %ZBL=zero baseline
corrected = zeros(size(data));

debug = zeros([size(data,1) size(data,2)]);

for i=1:size(data,1)
    for j=1:size(data,2)
        leftBL = mean(data(i,j,1:5));
        rightBL = mean(data(i,j,size(specData,1)-5:end));
        meanbaseline = mean([leftBL rightBL]);
        
        pixelrange = meanbaseline-min(data(i,j,:)); %difference between baseline and the minimum value in pixel spectrum
        
        if pixelrange == 0
            corrected(i,j,:) = nan;
        end
        
        debug(i,j) = pixelrange;
        
        for k=1:size(data,3)
            corrected(i,j,k) = (baselinerange / pixelrange) * (data(i,j,k) - meanbaseline) - globalFraction*specZBL(k,1) + globalmeanBL;
        end
    end
end
