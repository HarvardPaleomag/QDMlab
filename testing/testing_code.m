%expData = load('/Volumes/backups/antares/mike/MIL/MIL3/FOV1/1000G/final_fits_(2x2).mat');
expData = load('Z:\antares\mike\MIL\MIL3\FOV1\NRM\14G\4x4Binned\final_fits_(4x4).mat');
expData = load('Z:\antares\mike\MIL\MIL3\FOV1\NRM\4x4Binned\final_fits_(4x4).mat');
% expData = load('/Volumes/backups/antares/tiled_blanks/mariner/FOV_2_2/final_fits_(4x4).mat');
% expData = load('/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/blanks/2021_05_11/final_fits_(4x4).mat');
% expData = load('/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/speleothems/10-20-10/final_fits_(4x4).mat');
b = expData.B111ferro;
%%  
close all 
clc
[d1 d2 d3 d4] = B111(expData, 'fieldsPlot', 1);
%%
close all
clc
f = figure('units','normalized','outerposition',[0.15 0.4 0.6 0.6],'NumberTitle', 'off', 'Name', '');

ax1 = subplot(1,1,1);
[f,ax,im] = QDM_figure(expData.ledImg, 'ax', ax1, 'led', true, 'xc', 1:960, 'yc', 1:600, 'title', 'LED');
hold on

QDM_figure(b, 'alpha', 1e-13, 'ax', ax1, 'xc', 1:960, 'yc', 1:600, 'title', 'QDM')

% QDM_figure(b, 'alpha', 1, 'ax', ax)
%%
function [b111Ferro, b111Para, centerShift, dDelta] = B111(expData, plots)
    arguments
        expData
        plots.resonanceFieldPlot {mustBeBoolean(plots.resonanceFieldPlot)} = false
        plots.centerShiftPlot {mustBeBoolean(plots.centerShiftPlot)} = false
        plots.fieldsPlot {mustBeBoolean(plots.fieldsPlot)} = false
    end
    
    gamma = 0.0028;
    
    resPos1 = expData.leftPos.resonance;
    resPos2 = expData.rightPos.resonance;
    
    resNeg1 = expData.leftNeg.resonance;
    resNeg2 = expData.rightNeg.resonance;
    
    medianResPos = median((resPos2+resPos1)/2, 'all', 'omitnan');
    medianResNeg = median((resNeg2+resNeg1)/2, 'all', 'omitnan');
    med = mean([medianResPos medianResNeg]);
    
    %% shift away from the mean
    % negative fields
    dNeg1 = (med-resNeg1)/gamma;
    dNeg2 = (resNeg2-med)/gamma;
    % positive fields
    dPos1 = (med-resPos1)/gamma;
    dPos2 = (resPos2-med)/gamma;
    
    %% mean shift (i.e. mean(field) with pos/neg field applied
    meanPos =  (dPos1+dPos2)/2;
    meanNeg =  (dNeg1+dNeg2)/2;
    
    %% ferro / para component calculation
    b111Ferro = (meanPos - meanNeg)/2; 
    b111Para  = (meanPos + meanNeg)/2;
    
    %% center shift in pos/neg field direction
    deltaPos = (dPos2-dPos1)/2;
    deltaNeg = (dNeg2-dNeg1)/2;
    
    %% average center shift
    centerShift = (deltaPos + deltaNeg)/2;
    
    %% differences between pos/ neg field
    dDelta = (deltaNeg - deltaPos)/2;
    
    %% PLOTTING
    %% center shift
    if isequal(plots.centerShiftPlot, true)
        f = figure('units','normalized','outerposition',[0.15 0.4 0.7 0.4],'NumberTitle', 'off', 'Name', 'Center Shifts');
        ax1 = subplot(1,2,1);
        QDM_figure(deltaPos, 'ax', ax1, 'title', 'B^+');
        ax2 = subplot(1,2,2);
        QDM_figure(deltaNeg, 'ax', ax2, 'title', 'B^-');
        linkaxes([ax1 ax2]);
    end
    
    %% fields
    if isequal(plots.fieldsPlot, true)
        f = figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8],'NumberTitle', 'off', 'Name', 'B111ferro, B111para');
        ax1 = subplot(2,2,1);
        QDM_figure(b111Ferro, 'st', 10, 'preThreshold', 10, 'ax', ax1, 'title', 'B_{111} ferro');
        ax2 = subplot(2,2,2);
        QDM_figure(b111Para-median(b111Para, 'all'), 'st', 10, 'preThreshold', 10, 'ax', ax2, 'title', 'B_{111} para', 'mustBe', false);
        ax3 = subplot(2,2,3);
        QDM_figure(centerShift, 'st', 20, 'preThreshold', 1, 'ax', ax3, 'title', 'mean(center shift)');
        ax4 = subplot(2,2,4);
        QDM_figure(dDelta, 'st', 2, 'preThreshold', 1, 'ax', ax4, 'title', 'CS^- - CS^+');
        linkaxes([ax1 ax2 ax3 ax4]);
    end
    
    %% resonance fields
    if isequal(plots.resonanceFieldPlot, true)
        f = figure('units','normalized','outerposition',[0.15 0.4 0.7 0.4],'NumberTitle', 'off', 'Name', 'Resonance Fields');
        ax1 = subplot(2,2,1);
        QDM_figure(dPos1-median(dPos1,'all'), 'ax', ax1, 'title', 'B^+_1');
        ax2 = subplot(2,2,2);
        QDM_figure(dPos2-median(dPos2,'all'), 'ax', ax2, 'title', 'B^+_2');
        ax3 = subplot(2,2,3);
        QDM_figure(dNeg1-median(dNeg1,'all'), 'ax', ax3, 'title', 'B^-_1');
        ax4 = subplot(2,2,4);
        QDM_figure(dNeg2-median(dNeg2,'all'), 'ax', ax4, 'title', 'B^-_2');
        linkaxes([ax1 ax2 ax3 ax4]);
    end
end