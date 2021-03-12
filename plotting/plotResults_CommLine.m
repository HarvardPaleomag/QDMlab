function [] = plotResults_CommLine(dataFolder, type, kwargs)

arguments
    dataFolder
    type
    kwargs.checkPlot (1,1) {mustBeMember(kwargs.checkPlot, [1, 0])} = 0
end

close all;
 
myDir = dataFolder;

ledFiles = dir(fullfile(myDir,'led.csv'));   %grab the first / only CSV
ledImgPath = fullfile(myDir, ledFiles(1).name);

gamma = 0.0028;
 
if type == 'np  '
    negB111Output = load(fullfile(myDir, 'run_00000.matdeltaBFit.mat'));
    posB111Output = load(fullfile(myDir, 'run_00001.matdeltaBFit.mat'));
    
    negDiff = - real( (negB111Output.Resonance2-negB111Output.Resonance1)/2 / gamma );
    posDiff = real( (posB111Output.Resonance2-posB111Output.Resonance1)/2 / gamma );

    B111ferro = (posDiff + negDiff)/2;
    B111para = (posDiff - negDiff)/2;
else
  if type == 'nppn'
    negB111Output = load(fullfile(myDir, 'run_00000.matdeltaBFit.mat'));
    posB111Output = load(fullfile(myDir, 'run_00001.matdeltaBFit.mat'));
    
    negDiff = - real( (negB111Output.Resonance2-negB111Output.Resonance1)/2 / gamma );
    posDiff = real( (posB111Output.Resonance2-posB111Output.Resonance1)/2 / gamma );
    
    negB111Output2 = load(fullfile(myDir, 'run_00003.matdeltaBFit.mat'));
    posB111Output2 = load(fullfile(myDir, 'run_00002.matdeltaBFit.mat'));
    
    negDiffR = - real( (negB111Output2.Resonance2-negB111Output2.Resonance1)/2 / gamma );
    posDiffR = real( (posB111Output2.Resonance2-posB111Output2.Resonance1)/2 / gamma );
    
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

fitFailed = posB111Output.fitFailed | negB111Output.fitFailed;
rng = .03;

ledImg = load(ledImgPath);



%% SAVE results for plotting later
B111dataToPlot.negDiff = double(negDiff); B111dataToPlot.posDiff = double(posDiff); 
B111dataToPlot.B111ferro = double(B111ferro); B111dataToPlot.B111para = double(B111para);
B111dataToPlot.chi2Pos1 = double(chi2Pos1); B111dataToPlot.chi2Pos2 = double(chi2Pos2); 
B111dataToPlot.chi2Neg1 = double(chi2Neg1); B111dataToPlot.chi2Neg2 = double(chi2Neg2);
B111dataToPlot.ledImg = ledImg; B111dataToPlot.fitFailed = fitFailed; 
save(fullfile(myDir, 'B111dataToPlot.mat'), '-struct', 'B111dataToPlot');

%% PLOTS
r1 = nanmean(B111para(fitFailed))-rng;    r2 = nanmean(B111para(fitFailed))+rng;
 
f1=figure; imagesc( (negDiff) ); axis equal tight; caxis([-r2 -r1]); colorbar; colormap jet; title('Negative current B_{111} (gauss)'); set(gca,'YDir','normal');
f2=figure; imagesc( (posDiff) ); axis equal tight; caxis([r1 r2]); colorbar; colormap jet; title('Positive current B_{111} (gauss)'); set(gca,'YDir','normal');
f3=figure; imagesc( B111ferro ); axis equal tight; caxis(-.1 + [-rng rng]); colorbar; colormap jet; title('Positive + negative ferro B_{111} (gauss)'); set(gca,'YDir','normal');
f4=figure; imagesc(ledImg); axis equal tight; colorbar; colormap gray; caxis auto; title('LED image'); set(gca,'YDir','normal');
f5=figure; imagesc( B111para ); axis equal tight; caxis([r1 r2]); colorbar; colormap jet; title('Positive + negative para B_{111} (gauss)'); set(gca,'YDir','normal');
 
f6=figure('Name','data','units','normalized','outerposition',[0 0 1 1]); set(gca,'YDir','normal');
s1 = subplot(2,2,1); imagesc( (negDiff) ,'hittest', 'off'); axis equal tight; caxis auto; colorbar; colormap(s1,jet); title('Negative current B_{111} (gauss)'); set(gca,'YDir','normal');
s2 = subplot(2,2,2); imagesc( (posDiff) ,'hittest', 'off'); axis equal tight; caxis auto; colorbar; colormap(s2,jet); title('Positive current B_{111} (gauss)'); set(gca,'YDir','normal');
s3 = subplot(2,2,3); imagesc( B111ferro ,'hittest', 'off'); axis equal tight;  caxis(mean2(B111ferro) + [-rng rng]); colorbar; colormap(s3,jet); title('Positive + negative ferro B_{111} (gauss)'); set(gca,'YDir','normal');
s4 = subplot(2,2,4); imagesc( (ledImg) ,'hittest', 'off'); axis equal tight; colorbar; colormap(s4,gray); caxis auto; title('LED image'); set(gca,'YDir','normal');
sgtitle(', B111 points up and out of page');
ax = [s1,s2,s3];
linkaxes(ax);


if kwargs.checkPlot
    fig = figure('Name', 'spectra');
    bin.ledImg = ledImg;
    bin.B111ferro = B111ferro;
    binSize = detect_binning(bin);
    points = [0,0,0];
    for s = ax
       set(s,'ButtonDownFcn',{@clickFitFcn, binSize, ...
           negB111Output, negB111Output, points, ax, f6})
    end
    waitfor(f6)
    close all
    drawnow;
end

saveas(f1, fullfile(myDir, 'negCurrent.png'),'png');
saveas(f2, fullfile(myDir, 'posCurrent.png'),'png');
saveas(f3, fullfile(myDir, 'ferromagImg.png'),'png');
saveas(f4, fullfile(myDir, 'ledImg.png'),'png');
saveas(f5, fullfile(myDir, 'paramagImg.png'),'png');
saveas(f6, fullfile(myDir, 'allPlots.png'),'png');
end
function clickFitFcn(hObj, event, binSize, posB111Output, negB111Output, points, ax, fig)
    % Get click coordinate
    spec = findobj( 'Type', 'Figure', 'Name', 'spectra' );
    dat = findobj( 'Type', 'Figure', 'Name', 'data' );
    
    click = event.IntersectionPoint;
    x = round(click(1));
    y = round(click(2));

    titleTxt = sprintf('X: %4i (%4i) Y: %4i (%4i)', ...
        round(x),round(x)*binSize,round(y), round(y)*binSize);
    
    %% plot spectra
    set(0, 'currentfigure', spec)
    ax1 = subplot(1,2,1); cla(); hold on
    
    plot(ax1, posB111Output.Freqs1, squeeze(posB111Output.Resonance1(y,x,:)), 'k.-','DisplayName','+')
    plot(ax1, negB111Output.Freqs1, squeeze(negB111Output.Resonance1(y,x,:)), 'k.-','DisplayName','-')
    
    ax2 = subplot(1,2,2); cla(); hold on
    plot(ax2, posB111Output.Freqs2, squeeze(posB111Output.Resonance2(y,x,:)), 'k.-','DisplayName','+')
    plot(ax2, negB111Output.Freqs2, squeeze(negB111Output.Resonance2(y,x,:)), 'k.-','DisplayName','-')
    
    %% plot points
    set(0, 'currentfigure', dat)
    for a = [ax1 ax2]
        ylabel(a, 'Intensity')
        xlabel(a, 'f (Hz)')
        legend(a, 'Location','southwest', 'NumColumns',3)
    end
    
    for i = 1:3
        point = points(i);
        if point ~= 0
            delete(point);
        end
        set(dat, 'currentaxes', ax(i))
        hold on
        point = scatter(round(x),round(y),'xr');
        points(i) = point;
    end
    title(ax1, titleTxt)
end