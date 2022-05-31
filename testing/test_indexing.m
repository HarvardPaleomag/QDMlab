clear all
close all
clc
logLevel = 'debug';
expData = load("C:\Users\micha\OneDrive\Desktop\NRM\FOV1\run_00000.mat");
header = read_header("C:\Users\micha\OneDrive\Desktop\NRM\FOV1\run_00000_header.txt");
%%
[binDataNorm, freq] = prepare_raw_data(expData, header, 1, 1, normalize=false);
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
d = [];
for r= 1:size(A,1)
    for c= 1:size(A,2);
        i = xy2index(r,c, size(A),'type', 'binDataNorm');
        d = [d;[A(r,c), a(i)]];
    end
end
all(d(:,1)==d(:,2))
%% simpler test
clc
aux = a([1 2 3 4 5 6 7 8 9 10]);
d = [];
for i = 1:10
    [r,c] = index2xy(i, size(A),'type', 'binDataNorm');
    [i,r,c]
    d = [d A(r,c)];
end
aux
d
%%
clc
d = []
for i = 1:size(a,2)
    [r,c] = index2xy(i, size(A),'type', 'binDataNorm')
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
