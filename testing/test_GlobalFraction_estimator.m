expData = load(['C:\Users\micha\OneDrive\Desktop\NRM\FOV1\run_00000.mat']);
header = read_header('C:\Users\micha\OneDrive\Desktop\NRM\FOV1\run_00000_header.txt');
%%
globalFraction_estimator('expData', expData, 'header', header)
%%
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperPosition', [1 1 28 19]);
set(h,'PaperPositionMode','auto');
print(gcf, '-dpdf', '/Users/mike/Dropbox/science/_projects/QDM/figs/GlobalFraction_estimator.pdf','-bestfit');