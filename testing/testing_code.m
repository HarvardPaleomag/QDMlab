expData = load('/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/ALH84001/S10/FOV1_VISC_20-20/4x4Binned/final_fits_(4x4).mat');
% expData = load('/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/blanks/2021_05_11/final_fits_(4x4).mat');
% expData = load('/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/speleothems/10-20-10/final_fits_(4x4).mat');
b = expData.B111ferro;
%%
close all 
clc
[d1 d2 d3 d4] = B111(expData);

f = figure('units','normalized','outerposition',[0.15 0.4 0.7 0.4],'NumberTitle', 'off', 'Name', 'title');

ax1 = subplot(1,2,1);
QDM_figure(d3, 'ax', ax1, 'title', 'mean(delta_{res})', 'scalebar', 250);

ax2 = subplot(1,2,2);
QDM_figure(d4, 'ax', ax2, 'st', 2, 'title', 'delta_{field}');
linkaxes([ax1 ax2]);
% d1 = filter_hot_pixels(d1, 'threshold', 10);
% d2 = filter_hot_pixels(d2, 'threshold', 10);
% d3 = filter_hot_pixels(d3, 'threshold', 5);

f = figure('units','normalized','outerposition',[0.05 0.4 0.9 0.35],'NumberTitle', 'off', 'Name', 'This is the figure title');
ax1 = subplot(1,3,1);
QDM_figure(d1, 'st', 10, 'preThreshold', 10, 'ax', ax1, 'title', 'B_{111} ferro');
ax2 = subplot(1,3,2);
QDM_figure(d2-median(d2, 'all'), 'st', 10, 'preThreshold', 10, 'ax', ax2, 'title', 'B_{111} para', 'mustBe', false);
ax3 = subplot(1,3,3);
QDM_figure(d3, 'st', 10, 'preThreshold', 0.02, 'ax', ax3, 'title', 'delta B_{111}');
linkaxes([ax1 ax2 ax3]);


%%
TF = isoutlier(d1,'percentiles',[0,99.5]);
numel(nonzeros(TF))
% QDM_figure(d1.*abs(TF-1));
f = figure('units','normalized','outerposition',[0.2 0.4 0.6 0.3],'NumberTitle', 'off', 'Name', 'This is the figure title');
ax1 = subplot(1,2,1);
QDM_figure(d1,'ax', ax1);
ax2 = subplot(1,2,2);
QDM_figure(filloutliers(d1, 'linear', 'gesd','ThresholdFactor', 5e-14),'ax', ax2)
linkaxes([ax1 ax2]);
%%
function [b111Ferro, b111Para, delta, dDelta] = B111(expData)
    gamma = 0.0028;
    
    resPos1 = expData.leftPos.resonance;
    resPos2 = expData.rightPos.resonance;
    
    resNeg1 = expData.leftNeg.resonance;
    resNeg2 = expData.rightNeg.resonance;
    
    medianResPos = median((resPos2+resPos1)/2, 'all', 'omitnan');
    medianResNeg = median((resNeg2+resNeg1)/2, 'all', 'omitnan');
    med = mean([medianResPos medianResNeg]);
    
    dPos1 = (med-resPos1)/gamma;
    dPos2 = (resPos2-med)/gamma;
    dNeg1 = (med-resNeg1)/gamma;
    dNeg2 = (resNeg2-med)/gamma;
    
    meanPos =  (dPos1+dPos2)/2;
    meanNeg =  (dNeg1+dNeg2)/2;
    
    b111Ferro = (meanPos - meanNeg)/2; 
    b111Para  = (meanPos + meanNeg)/2;
    
    deltaPos = (dPos2-dPos1)/2;
    deltaNeg = (dNeg2-dNeg1)/2;
    
    delta = (deltaPos + deltaNeg)/2;
    dDelta = (deltaNeg - deltaPos)/2;

%     b111Ferro = (posDiff + negDiff)/2;
    f = figure('units','normalized','outerposition',[0.15 0.4 0.7 0.4],'NumberTitle', 'off', 'Name', 'title');

    ax1 = subplot(2,2,1);
    QDM_figure(dPos1-median(dPos1,'all'), 'ax', ax1, 'title', 'dPos1');

    ax2 = subplot(2,2,2);
    QDM_figure(dPos2-median(dPos2,'all'), 'ax', ax2, 'title', 'dPos2');
    
    ax3 = subplot(2,2,3);
    QDM_figure(dNeg1-median(dNeg1,'all'), 'ax', ax3, 'title', 'dNeg1');

    ax4 = subplot(2,2,4);
    QDM_figure(dNeg2-median(dNeg2,'all'), 'ax', ax4, 'title', 'dNeg1');
    linkaxes([ax1 ax2 ax3 ax4]);
end

% negDiff = - real( (negRes.Resonance2-negRes.Resonance1)/2 / gamma );
% posDiff =   real( (posRes.Resonance2-posRes.Resonance1)/2 / gamma );
% 
% B111ferro = (posDiff + negDiff)/2;
% B111para = (posDiff - negDiff)/2;