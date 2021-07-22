clear all
clc
fits = load('/Volumes/exchange/mike/composite_test/polymer_spherules_Mount2_FOV1_IRM10000G/spherule_allGF_2.mat');
%%
close all
clc

GFs = fits.GFs;
nElem = 500;
GF = zeros(1,nElem);

yStop = size(fits.resonance(:,:,1),1)/3;
xStop = size(fits.resonance(:,:,1),2)/3;

slope = 0.0106; %determined from fitting Chi2 values of blank map
% zfs = 2.870;
res = fits.resonance(:,:,1);
zfs = median(res, 'all');
res(1:yStop,:) = zfs;
res(2*yStop:end,:) = zfs;
res(:,1:xStop) = zfs;
res(:,2*xStop:end) = zfs;
% % imagesc(res)

GF0 = reshape(res, 1,[]);
chi0 = reshape(fits.chiSquares(:,:,1),1,[]);
[vals, idx] = sort(abs(zfs-GF0),'descend');
idx = idx(1:nElem);

for i = 2:size(GFs,2)
    chi = reshape(fits.chiSquares(:,:,i), 1,[]);
    chi = chi(idx);
    chiSmaller = chi - (slope * GFs(i)) < chi0(idx);

    GF(chiSmaller) = GFs(i);
end
% globalFraction = mean(GF(GF > 0));
globalFraction = mean(GF);

[f,ax] = QDM_figure(fits.resonance(:,:,1));
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
aux = res;
aux(res == zfs) = 0;
visboundaries(aux,'Color','b')
figure; histogram(GF)
%%
globalFraction_estimator('expData','/Volumes/backups/antares/roger/diagnostics/polymer_spherules/Mount2_FOV1_IRM2/IRM10000G');
%%
[globalFraction, GF, idx] = determine_globalFraction(data)