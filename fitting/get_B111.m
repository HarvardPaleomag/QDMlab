function fit = get_B111(dataFolder, folderName, type, fit, binSize, kwargs)
%[fits] = plotResults_CommLine(dataFolder, folderName, type, fits, binSize; 'checkPlot')

arguments
    dataFolder
    folderName
    type
    fit
    binSize
    kwargs.checkPlot (1, 1) {mustBeBoolean(kwargs.checkPlot)} = false;
    kwargs.crop = 'none'
end

close all;

ledFiles = dir(fullfile(dataFolder,'led.csv'));   %grab the first / only CSV
ledImgPath = fullfile(dataFolder, ledFiles(1).name);
ledImg = load(ledImgPath);

laserFiles = dir(fullfile(dataFolder,'laser.csv'));   %grab the first / only CSV
laserImgPath = fullfile(dataFolder, laserFiles(1).name);
laserImg = load(laserImgPath);

if ~strcmp(kwargs.crop, 'none')
    x0 = kwargs.crop(1);
    y0 = kwargs.crop(2);
    x1 = kwargs.crop(3)+x0;
    y1 = kwargs.crop(4)+y0;
    ledImg = ledImg(y0:y1,x0:y1);
    laserImg = laserImg(y0:y1,x0:x1);
    fit.crop = kwargs.crop
end
gamma = 0.0028;
 
if type == 'np  '
    negB111Output = load(fullfile(dataFolder, folderName, 'run_00000.matdeltaBFit.mat'));
    posB111Output = load(fullfile(dataFolder, folderName, 'run_00001.matdeltaBFit.mat'));
    
    negDiff = - real( (negB111Output.Resonance2-negB111Output.Resonance1)/2 / gamma );
    posDiff =   real( (posB111Output.Resonance2-posB111Output.Resonance1)/2 / gamma );

    B111ferro = (posDiff + negDiff)/2;
    B111para  = (posDiff - negDiff)/2;
else
  if type == 'nppn'
    negB111Output = load(fullfile(dataFolder, folderName, 'run_00000.matdeltaBFit.mat'));
    posB111Output = load(fullfile(dataFolder, folderName, 'run_00001.matdeltaBFit.mat'));
    
    negDiff = - real( (negB111Output.Resonance2-negB111Output.Resonance1)/2 / gamma );
    posDiff =   real( (posB111Output.Resonance2-posB111Output.Resonance1)/2 / gamma );
    
    negB111Output2 = load(fullfile(dataFolder, folderName, 'run_00003.matdeltaBFit.mat'));
    posB111Output2 = load(fullfile(dataFolder, folderName, 'run_00002.matdeltaBFit.mat'));
    
    negDiffR = - real( (negB111Output2.Resonance2-negB111Output2.Resonance1)/2 / gamma );
    posDiffR =   real( (posB111Output2.Resonance2-posB111Output2.Resonance1)/2 / gamma );
    
    B111ferro = (posDiff + negDiff + posDiffR + negDiffR)/4;  %must divide ferro part by 2
    B111para = (posDiff - negDiff + posDiffR - negDiffR)/4;
  end
end
 
% get Chi squared values for pos and neg. fields / left right fit
chi2Pos1 = posB111Output.chiSquares1;
chi2Pos2 = posB111Output.chiSquares2;
chi2Neg1 = negB111Output.chiSquares1;
chi2Neg2 = negB111Output.chiSquares2;
% reshape the chi2 to conform to pixels
chi2Pos1 = reshape(chi2Pos1, size(B111ferro));
chi2Pos2 = reshape(chi2Pos2, size(B111ferro));
chi2Neg1 = reshape(chi2Neg1, size(B111ferro));
chi2Neg2 = reshape(chi2Neg2, size(B111ferro));

%% SAVE results for plotting later
B111dataToPlot.negDiff = double(negDiff); 
B111dataToPlot.posDiff = double(posDiff); 

B111dataToPlot.B111ferro = double(B111ferro); 
B111dataToPlot.B111para = double(B111para);

B111dataToPlot.chi2Pos1 = double(chi2Pos1); 
B111dataToPlot.chi2Pos2 = double(chi2Pos2); 
B111dataToPlot.chi2Neg1 = double(chi2Neg1); 
B111dataToPlot.chi2Neg2 = double(chi2Neg2);

B111dataToPlot.ledImg = ledImg; 
B111dataToPlot.laser = laserImg;

%% add to final_fits structure
fit.negDiff = negDiff; fit.posDiff = posDiff; 
fit.B111ferro = B111ferro; fit.B111para = B111para;
fit.ledImg = ledImg; fit.laserImg = laserImg;

%% determine overall fit Success
pixelAlerts = posB111Output.pixelAlerts | negB111Output.pixelAlerts;
B111dataToPlot.pixelAlerts = pixelAlerts; 
fit.pixelAlerts = pixelAlerts;

fit.B111 = B111(fit);

save(fullfile(dataFolder, folderName, 'B111dataToPlot.mat'), '-struct', 'B111dataToPlot');
end