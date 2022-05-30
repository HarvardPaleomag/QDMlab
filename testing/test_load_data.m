% testing load_data
clear
logLevel = 'debug';
clc

% loading specified folder
folder = "C:\Users\micha\OneDrive\Desktop\NRM\FOV1";

%% only run_00000
f1 = load_data(folder, 1);

%% only run_00001
f2 = load_data(folder, 2);

%% all run files
f3 = load_data(folder);

%% loading data without any specification
% should open dialog to select a folder with run_** files
no_spec = load_data();