%% B111 test
expData = load('/Users/mike/Desktop/ALH84001/S10/FOV1_VISC_20-20/4x4Binned/B111dataToPlot.mat');
%%
% should return true
is_B111(expData) == true
%%
% should return true, and the names
[bool, dataName, ledName] = is_B111(expData);
1 == bool && all(['B111ferro' 'ledImg'] == [dataName, ledName])

%% Bz test
expData = load('/Users/mike/Desktop/ALH84001/S10/FOV1_VISC_20-20/4x4Binned/Bz_uc0.mat');
%%
% should return true
is_B111(expData) == false
%%
% should return true, and the names
[bool, dataName, ledName] = is_B111(expData);
0 == bool && all(['Bz' 'newLED'] == [dataName, ledName])