clear all
close all
clc
logLevel = 'debug';
expData = load('/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/QDM_data/NRM_MIL/run_00000.mat');
%%
[binDataNorm, freq] = prepare_raw_data(expData, 1, 1);
%%
data = QDMreshape(expData.imgStack1, expData.imgNumRows, expData.imgNumCols);
%%
clc
all(binDataNorm-data < 0.01, 'all') % only if normalization isnt done in prepare_raw_data
%%
a = 1:50;
A = QDMreshape(a, 5,10);
gpudata = reshape(A, size(a,2), 1); % make it into 2d matrix
% gpudata = transpose(A); %transpose to make it 51 x pixels
%% test binDataNorm arrays
d = []
for r= 1:size(A,1)
    for c= 1:size(A,2);
        i = xy2index(r,c, size(A),'type', 'binDataNorm');
        d = [d,[a(i); A(r,c)]];
    end
end
d'
%%
clc
d = []
for i = 1:size(a,2)
    [r,c] = index2xy(i, size(A));
    d = [d,[a(i); A(r,c)]];
end
all(d(1,:) == d(2,:))

%% test gpudata arrays
d = []
for r= 1:size(A,1)
    for c= 1:size(A,2);
        i = xy2index(r,c, size(A), 'type', 'gpu');
        d = [d,[gpudata(i); A(r,c)]];
    end
end
all(d(1,:) == d(2,:))
%%
clc
d = []
for i = 1:size(a,2)
    [r,c] = index2xy(i, size(A), 'type', 'gpu');
    d = [d,[gpudata(i); A(r,c)]];
end
all(d(1,:) == d(2,:))
