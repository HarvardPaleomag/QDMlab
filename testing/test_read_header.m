% testing read_header function
clear
logLevel = 'debug';
clc

% reading only a file name
headerFile = "C:\Users\micha\OneDrive\Desktop\NRM\FOV1\run_00000_header.txt";
h1 = read_header(headerFile);

% reading a dir(file) structure
headerFile = dir(headerFile);
h2 = read_header(headerFile);

%reading without specifying path
h3 = read_header();

