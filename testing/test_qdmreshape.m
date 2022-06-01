%% testing QDMreshape
clear all
close all
clc
logLevel = 'debug';
expData = load("C:\Users\micha\OneDrive\Desktop\NRM\FOV1\run_00000.mat");
header = read_header("C:\Users\micha\OneDrive\Desktop\NRM\FOV1\run_00000_header.txt");
%%
[binDataNorm, freq] = prepare_raw_data(expData, 1, 1, header, normalize=false);
%%
data = QDMreshape(expData.imgStack1, expData.imgNumRows, expData.imgNumCols);
%% test if QDMreshape correctly reshapes the data
clc
all(binDataNorm-data < 0.01, 'all') % only if normalization isnt done in prepare_raw_data
%%
dataReshape = QDMreshape_reverse(data, expData.numFreqs);
all(dataReshape-expData.imgStack1 < 0.01, 'all') % only if normalization isnt done in prepare_raw_data
