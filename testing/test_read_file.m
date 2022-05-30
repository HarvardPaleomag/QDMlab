% testing file read_data
clear
logLevel = 'debug';
clc

% reading only a file name
file = "C:\Users\micha\OneDrive\Desktop\NRM\FOV1\run_00000.mat";
a = read_file(file);

% reading a dir(file) structure
runFile = dir(file);
b = read_file(runFile);

%reading without specifying path
c = read_file();