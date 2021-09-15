clear
logLevel = 'debug';
folders = {'/Users/mike/Dropbox/science/harvard/QDM2_data/antares/mike/MIL/MIL3/FOV1/NRM/2x2Binned';
           '/Users/mike/Dropbox/science/harvard/QDM2_data/antares/mike/MIL/MIL3/FOV1/150G/2x2Binned';
           '/Users/mike/Dropbox/science/harvard/QDM2_data/antares/mike/MIL/MIL3/FOV1/280G_2/2x2Binned';
           '/Users/mike/Dropbox/science/harvard/QDM2_data/antares/mike/MIL/MIL3/FOV1/400G/2x2Binned';
           '/Users/mike/Dropbox/science/harvard/QDM2_data/antares/mike/MIL/MIL3/FOV1/500G/2x2Binned';
           '/Users/mike/Dropbox/science/harvard/QDM2_data/antares/mike/MIL/MIL3/FOV1/600G/2x2Binned';
           '/Users/mike/Dropbox/science/harvard/QDM2_data/antares/mike/MIL/MIL3/FOV1/1000G/2x2Binned'};  
fname = 'B111dataToPlot.mat';
%%
res = demag_behavior(folders, 'fileName', fname, 'bootStrapN', 1000, 'pixelShift', 3, 'filterProps', struct('threshold', 1.5));%, 'filterProps', struct('threshold', 5, 'chi', true, 'cutOff', 14))
%%
demag_behavior_plot(res, 'steps', [0,150,280,400,500,600,1000],mean=true)