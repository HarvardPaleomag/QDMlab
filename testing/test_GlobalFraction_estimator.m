expData = load('/Users/mike/Desktop/NRM/run_00000.mat');
%%
globalFraction_estimator(expData)
%%
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperPosition', [1 1 28 19]);
set(h,'PaperPositionMode','auto');
print(gcf, '-dpdf', '/Users/mike/Dropbox/science/_projects/QDM/figs/GlobalFraction_estimator.pdf','-bestfit');