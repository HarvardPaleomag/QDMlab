clear all
clc
fits = load('/Volumes/exchange/mike/composite_test/polymer_spherules_Mount2_FOV1_IRM10000G/spherule_allGF_2.mat');
%%
close all
clc

GFs = fits.GFs;
nElem = 500;
nElem = numel(fits.resonance(:,:,1))/800;
% nElem = numel(fits.resonance(:,:,1));
% fits.resonance(fits.resonance>2.8446) = nan;
GF = zeros(1,nElem);

yStop = 100;size(fits.resonance(:,:,1),1)/3;
xStop = 200;size(fits.resonance(:,:,1),2)/3;

slope = 0.0176; %determined from fitting Chi2 values of blank map
% zfs = 2.870;
res = fits.resonance(:,:,1);
zfs = median(res, 'all', 'omitnan');
% res(1:yStop,:) = zfs;
% res(2*yStop:end,:) = zfs;
% res(:,1:xStop) = zfs;
% res(:,2*xStop:end) = zfs;
% % imagesc(res)

sortVals = zeros(size(res));
% sortVals(yStop:2*yStop, xStop:2*xStop) = fits.chiSquares(yStop:2*yStop, xStop:2*xStop,1);
% sortVals = fits.chiSquares(:,:,1);
sortVals(yStop:end-yStop, xStop:end-xStop) = abs(zfs-fits.resonance(yStop:end-yStop, xStop:end-xStop,1));
% sortVals(yStop:end-yStop, xStop:end-xStop) = abs(zfs-fits.contrastC(yStop:end-yStop, xStop:end-xStop,1));

% sortVals(yStop:2*yStop, xStop:2*xStop) = fits.resonance(yStop:2*yStop, xStop:2*xStop,1);
% sortVals(yStop:2*yStop, xStop:2*xStop) = abs(sortVals(yStop:2*yStop, xStop:2*xStop));
% sortVals = fits.resonance(:,:,1);
sortVals(fits.contrastA(:,:,1) > 0.9) = 0;
sortVals(fits.contrastB(:,:,1) > 0.9) = 0;
sortVals(fits.contrastC(:,:,1) > 0.9) = 0;

GF0 = reshape(res, 1,[]);

chi0 = reshape(fits.chiSquares(:,:,1),1,[]);

% [vals, idx] = sort(abs(zfs-GF0),'descend');
[vals, idx] = sort(reshape(sortVals, 1, []),'descend');

idx = idx(1:nElem);

pquad = [0.0350, 0.0001, 0.0216]; % quadratic fit for blank map

p5 = [0.0176, 0.0205]; % linear fit for blank map (-GF 0.5)
p6 = [0.0243, 0.0194]; % linear fit for blank map (-GF 0.6)

slope = p5(1);

for i = 2:size(GFs,2)
    chi = reshape(fits.chiSquares(:,:,i), 1,[]);
    chi = chi(idx);
    chiSmaller = chi - (slope * GFs(i)) < chi0(idx);
%     chiSmaller = chi - (p(1) * GFs(i)^2 + p(2) * GFs(i)) < chi0(idx);
    GF(chiSmaller) = GFs(i);
end
globalFraction = mean(GF(GF > 0));
% globalFraction = mean(GF);

[f,ax] = QDM_figure(fits.resonance(:,:,1), std=2);
hold on

d = fits.resonance(:,:,1);
for i = 1:numel(d)
    if ~ismember(idx,i)
        d(i) = 0;
    else
        d(i) = 10;
    end
end
xc = 1:size(fits.chiSquares(:,:,1),2);
yc = 1:size(fits.chiSquares(:,:,1),1);

imAlpha=ones(size(d));
imAlpha(isnan(d)) = 0;
im = visboundaries(ax, d);
aux = sortVals;
% aux(res == zfs) = 0;
visboundaries(aux,'Color','b')
figure; histogram(GF, binedges=fits.GFs)
xlabel('GF')
ylabel('n')
%%
close all
GFs = fits.GFs;
f1 = polyval(pquad,fits.GFs)';
f1 = f1 - min(f1);
f1 = f1 / max(f1);
% meanChi = squeeze(mean(fits.chiSquares, [1,2]));
meanChi = squeeze(mean(fits.chiSquares(yStop:2*yStop, xStop:2*xStop,:), [1,2]));
meanChi = meanChi - min(meanChi);
meanChi = meanChi / max(meanChi);

plot(GFs, f1)
hold on
plot(GFs, meanChi)
plot(GFs, meanChi-f1)
%%
close all
i = 1;
h = figure;
% axis tight manual % this ensures that getframe() returns a consistent size
filename = 'animated_contrasts.gif';
for i = 1:size(GFs,2)
    cla
    hold on
%     plot(reshape(fits.width(:,:,i),[],1),reshape(fits.chiSquares(:,:,i),[],1),  '.')
    plot(reshape(fits.contrastA(:,:,i),[],1),reshape(fits.chiSquares(:,:,i),[],1),  '.')
    plot(reshape(fits.contrastB(:,:,i),[],1),reshape(fits.chiSquares(:,:,i),[],1),  '.')
    plot(reshape(fits.contrastC(:,:,i),[],1),reshape(fits.chiSquares(:,:,i),[],1),  '.')
%     [GFs(i), numel(nonzeros(fits.width(:,:,i)>4.1093e-04*2))]
%     [GFs(i), numel(nonzeros(fits.contrastC(:,:,i)<0))]
%     xlim([0,0.1])
    legend('A','B','C')
    xlabel('contrast')
    ylabel('Chi^2')
    pause(0.2)
    drawnow
%     % Capture the plot as an image 
%     frame = getframe(h); 
%     im = frame2im(frame); 
%     [imind,cm] = rgb2ind(im,256); 
%     % Write to the GIF File 
%     if i == 1 
%       imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%     else 
%       imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%     end 
end
%%


a1 = subplot(2,2,1);
QDM_figure(fits.contrastA(:,:,1), 'std', 2, 'mustBe', 'pos', 'ax', a1)
a1 = subplot(2,2,2);
QDM_figure(fits.contrastB(:,:,1), 'std', 2, 'mustBe', 'pos', 'ax', a1)
a1 = subplot(2,2,3);
QDM_figure(fits.contrastC(:,:,1), 'std', 2, 'mustBe', 'pos', 'ax', a1)
    
    
%%
globalFraction_estimator('expData','/Volumes/backups/antares/roger/diagnostics/polymer_spherules/Mount2_FOV1_IRM2/IRM10000G');
%%
[globalFraction, GF, idx] = determine_globalFraction(data)