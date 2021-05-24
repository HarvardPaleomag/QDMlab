clear
logLevel = 'debug';
folders = {'/Users/mike/Dropbox/science/harvard/QDM2_data/antares/mike/MIL/MIL3/FOV1/NRM/2x2Binned',
           '/Users/mike/Dropbox/science/harvard/QDM2_data/antares/mike/MIL/MIL3/FOV1/150G/2x2Binned'
           '/Users/mike/Dropbox/science/harvard/QDM2_data/antares/mike/MIL/MIL3/FOV1/280G_2/2x2Binned'
           '/Users/mike/Dropbox/science/harvard/QDM2_data/antares/mike/MIL/MIL3/FOV1/400G/2x2Binned'
           '/Users/mike/Dropbox/science/harvard/QDM2_data/antares/mike/MIL/MIL3/FOV1/500G/2x2Binned'
           '/Users/mike/Dropbox/science/harvard/QDM2_data/antares/mike/MIL/MIL3/FOV1/600G/2x2Binned'
           '/Users/mike/Dropbox/science/harvard/QDM2_data/antares/mike/MIL/MIL3/FOV1/1000G/2x2Binned'};  
fname = 'B111dataToPlot.mat';
%%
estimate_coercivity(folders, 'fileName', fname, 'filterProps', struct('threshold', 1, 'chi', true))