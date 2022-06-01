clear all
close all
clc
logLevel = 'debug';
dataFile = "C:\Users\micha\OneDrive\Desktop\NRM\FOV1\run_00000.mat";
headerFile = "C:\Users\micha\OneDrive\Desktop\NRM\FOV1\run_00000_header.txt";

%% load data/header
expData = load(dataFile);
header = read_header(headerFile);
%%
[binDataNorm, freq] = prepare_raw_data(expData, 1, 1, header, normalize=false);
%%
data = QDMreshape(expData.imgStack1, expData.imgNumRows, expData.imgNumCols);
%%
gpudata = gpu_data_reshape(binDataNorm);
%% visualizing the different arrays
% if the indexing is correct all these plots should shopw the spectrum of
% the same pixel
close all
idx = 200;
hold on
% plot original spectrum of the unflattened array
plot(expData.imgStack1(:,idx))
% get index of the reshaped array
[row,col] = index2xy(idx, size(expData.imgStack1), type='binDataNorm');
gpuidx = xy2index(row, col, size(binDataNorm), 'gpu');
plot(squeeze(binDataNorm(row,col,:)),'x')
plot(squeeze(gpudata(:, gpuidx)), 'o')
legend('binDataNorm','index2xy','xy2index')


%% SIMPLE indexing testing
a = 1:50;
A = QDMreshape(a, 5,10);
gpudata = reshape(A, size(a,2), 1); % make it into 2d matrix
% gpudata = transpose(A); %transpose to make it 51 x pixels

%% test binDataNorm arrays
d = [];
for r= 1:size(A,1)
    for c= 1:size(A,2);
        i = xy2index(r,c, size(A),'type', 'raw');
        d = [d;[A(r,c), a(i)]];
    end
end
all(d(:,1)==d(:,2))
%% simpler test
clc
aux = a([1 2 3 4 5 6 7 8 9 10]);
d = [];
for i = 1:10
    [r,c] = index2xy(i, size(A),'type', 'raw');
    [i,r,c]
    d = [d A(r,c)];
end
aux
d
%%
clc
d = []
for i = 1:size(a,2)
    [r,c] = index2xy(i, size(A),'type', 'raw')
    d = [d;[a(i), A(r,c)]];
end
all(d(:,1) == d(:,2))

%% test gpudata arrays
clc
d = []
for r= 1:size(A,1)
    for c= 1:size(A,2);
        i = xy2index(r,c, size(A), 'type', 'gpu');
        d = [d;[A(r,c), gpudata(i)]];
    end
end
d
all(d(:,1) == d(:,2))
%%
clc
d = [];
for i = 1:size(a,2)
    [r,c] = index2xy(i, size(A), 'type', 'gpu');
    [i,r,c]
    d = [d;[gpudata(i), A(r,c)]];
end
d
all(d(:,1) == d(:,2))
