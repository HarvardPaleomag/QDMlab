function [] = plotResults_CommLine(INFILE,type)
 
close all;
 
myDir = INFILE;

ledFiles = dir(fullfile(myDir,'led.csv'));   %grab the first / only CSV
ledImgPath = fullfile(myDir, ledFiles(1).name);

gamma = 0.0028;
 
if type == 'np  '
  negB111Output = load(fullfile(myDir, 'run_00000.matdeltaBFit.mat'));
  posB111Output = load(fullfile(myDir, 'run_00001.matdeltaBFit.mat'));
  
  negDiff = - real( (negB111Output.Resonance2-negB111Output.Resonance1)/2 / gamma );
  posDiff = real( (posB111Output.Resonance2-posB111Output.Resonance1)/2 / gamma );
  
  B111ferro = (posDiff + negDiff)/2;  %must divide ferro part by 2
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
chi2Pos1 = posB111Output.chi_squares1;
chi2Pos2 = posB111Output.chi_squares2;
chi2Neg1 = negB111Output.chi_squares1;
chi2Neg2 = negB111Output.chi_squares2;
% reshape the chi2 to conform to pixels
chi2Pos1 = reshape(chi2Pos1, size(B111ferro));
chi2Pos2 = reshape(chi2Pos2, size(B111ferro));
chi2Neg1 = reshape(chi2Neg1, size(B111ferro));
chi2Neg2 = reshape(chi2Neg2, size(B111ferro));
 
rng = .03;
% ledImg = load( [ myDir 'LED_afterRun.csv'] );
ledImg = load(ledImgPath);
 
FitCvgBoth = negB111Output.FitCvg .* posB111Output.FitCvg;
 
r1 = mean2(B111para(FitCvgBoth==1))-rng;    r2 = mean2(B111para(FitCvgBoth==1))+rng;
 
f1=figure; imagesc( (negDiff) ); axis equal tight; caxis([-r2 -r1]); colorbar; colormap jet; title('Negative current B_{111} (gauss)'); set(gca,'YDir','normal');
f2=figure; imagesc( (posDiff) ); axis equal tight; caxis([r1 r2]); colorbar; colormap jet; title('Positive current B_{111} (gauss)'); set(gca,'YDir','normal');
% f3=figure; imagesc( B111ferro ); axis equal tight; caxis(mean2(B111ferro) + [-.5 .5]); colorbar; title('Positive + negative ferro part non-balanced B_{111} (gauss)');
f3=figure; imagesc( B111ferro ); axis equal tight; caxis(-.1 + [-rng rng]); colorbar; colormap jet; title('Positive + negative ferro B_{111} (gauss)'); set(gca,'YDir','normal');
f5=figure; imagesc( B111para ); axis equal tight; caxis([r1 r2]); colorbar; colormap jet; title('Positive + negative para B_{111} (gauss)'); set(gca,'YDir','normal');
 
% ledImg = flipud(ledImg);
f4=figure; imagesc(ledImg); axis equal tight; colorbar; colormap gray; caxis auto; title('LED image'); set(gca,'YDir','normal');
 
f6=figure('units','normalized','outerposition',[0 0 1 1]); set(gca,'YDir','normal');
s1 = subplot(2,2,1); imagesc( (negDiff) ); axis equal tight; caxis auto; colorbar; colormap(s1,jet); title('Negative current B_{111} (gauss)'); set(gca,'YDir','normal');
s2 = subplot(2,2,2); imagesc( (posDiff) ); axis equal tight; caxis auto; colorbar; colormap(s2,jet); title('Positive current B_{111} (gauss)'); set(gca,'YDir','normal');
s3 = subplot(2,2,3); imagesc( B111ferro ); axis equal tight;  caxis(mean2(B111ferro) + [-rng rng]); colorbar; colormap(s3,jet); title('Positive + negative ferro B_{111} (gauss)'); set(gca,'YDir','normal');
% s4 = subplot(2,2,4); imagesc( B111para ); axis equal tight; caxis([r1 r2]); colorbar; colormap jet; title('Positive + negative para B_{111} (gauss)'); set(gca,'YDir','normal');
s4 = subplot(2,2,4); imagesc( (ledImg) ); axis equal tight; colorbar; colormap(s4,gray); caxis auto; title('LED image'); set(gca,'YDir','normal');
sgtitle(', B111 points up and out of page');
 
 
 
%%
 
%%
%save(fullfile(myDir, 'B111dataToPlot.mat'],'negDiff','posDiff', 'B111ferro', 'B111para', 'ledImg', 'FitCvgBoth', 'chi2Pos1', 'chi2Pos2', 'chi2Neg1', 'chi2Neg2');
saveas(f1, fullfile(myDir, 'negCurrent.png'),'png');
saveas(f2, fullfile(myDir, 'posCurrent.png'),'png');
saveas(f3, fullfile(myDir, 'ferromagImg.png'),'png');
saveas(f4, fullfile(myDir, 'ledImg.png'),'png');
saveas(f5, fullfile(myDir, 'paramagImg.png'),'png');
saveas(f6, fullfile(myDir, 'allPlots.png'),'png');
 
%%
save(fullfile(myDir, 'B111dataToPlot.mat'),'negDiff','posDiff', 'B111ferro', 'B111para', 'ledImg', 'FitCvgBoth', 'chi2Pos1', 'chi2Pos2', 'chi2Neg1', 'chi2Neg2');
 
 
%%
% aaa=csvread('K:\Data\Analysis\20170210a_abom_GRA95229_B111_FOV2\afterRun2_LED_fullFOV_img0001.csv');
% figure; imagesc( aaa ); axis equal tight; title('LED image FOV1, big area'); colorbar; colormap gray; caxis([1000/4 7000/4]);  set(gca,'YDir','normal');
% hold on; plot([600 1700 600 1700], [600 600 1900 1900], 'r+');
