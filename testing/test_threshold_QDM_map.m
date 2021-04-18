d = '/Users/mike/Desktop/NRM/4x4Binned_GF/B111dataToPlot.mat';
%%
threshold_QDM_map(d, 'chi', true, 'cutOff', 5);
load('/Users/mike/Desktop/NRM/4x4Binned_GF/B111dataToPlot_thresh(5).mat');
QDM_figure(B111ferro);
%%
