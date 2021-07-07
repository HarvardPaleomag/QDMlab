clear all
close all
clc
logLevel = 'debug';
%%
nrm = load('/Users/mike/Dropbox/science/harvard/QDM2_data/mike/MIL/MIL3/2x2/FOV1/0NRM/B111dataToPlot.mat');
G400 = load('/Users/mike/Dropbox/science/harvard/QDM2_data/mike/MIL/MIL3/2x2/FOV1/600G/B111dataToPlot.mat');
%%
close all
f = figure('units','normalized','outerposition',[0.15 0.4 0.7 0.4],'NumberTitle', 'off', 'Name', 'title');

ax1 = subplot(1,2,1);
QDM_figure(nrm.ledImg, 'ax', ax1, 'led', true,'title','NRM');

ax2 = subplot(1,2,2);
QDM_figure(G400.ledImg, 'ax', ax2, 'led', true,'title', '600G');
linkaxes([ax1 ax2]);
%%
[tf, ref] = get_image_tform(nrm.ledImg, G400.ledImg, 'checkPlot', true);

%%
sem = imread('/Users/mike/Dropbox/science/_projects/miller_range_03346/data/102319_SEM_Mike/Electron Image SE/Electron Image SE_013.png');
sem = flip(sem);
sem = rot90(sem,1);
%%
subplot(1,2,1)
imagesc(sem)
a = subplot(1,2,2);
QDM_figure(nrm.ledImg, 'led', true, 'ax', a)
%%
[semTform semRef] = get_image_tform_complex(nrm.ledImg, sem, 'checkPlot',true);

%%
close all
clc
f = figure('units','normalized','outerposition',[0.15 0.4 0.7 0.4],'NumberTitle', 'off', 'Name', 'title');

semTransformed = tform_data(sem, semTform, semRef, 'binning', false);
a1 = subplot(1,2,1)
imagesc(a1, semTransformed)
axis equal xy
a2 = subplot(1,2,2)
QDM_figure(nrm.ledImg, 'led', true, 'ax', a2)
ylim(a1, a2.YLim)
set(a1, 'Position', [a1.Position(1:2) a2.Position(3:4)])

%%
close all
imshowpair(nrm.ledImg, semTransformed, 'Scaling', 'joint')