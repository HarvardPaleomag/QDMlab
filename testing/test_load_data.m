% testing load_data
clear
logLevel = 'debug';
clc

% loading data without any specification
% should open dialog to select a folder with run_** files
no_spec = load_data();