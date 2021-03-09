testData = load('/Users/mike/Desktop/NRM/run_00000.mat');
binSize = 2;
nRes = 1;
%%
[binDataNorm, freq] = prepare_raw_data(testData, binSize, nRes);
%%
corrected = correct_global(binDataNorm, 0.5);

%%
binSize=2;
i = 153; j = 134;
i_ = i*binSize-1:i*binSize-1+binSize-1;
j_ = j*binSize-1:j*binSize-1+binSize-1;

%% direct copy of old code
BackupExp = testData;

SpanX = 1:BackupExp.imgNumCols;
SpanY = 1:BackupExp.imgNumRows;

SpanXTrans = SpanX;
SpanYTrans = SpanY;

%Data loading
Fres1 = BackupExp.freqList(1:BackupExp.numFreqs)/1E9;   %everything is in GHz

Ares1stack = BackupExp.imgStack1;

%disp(peakwidth);
sweeplength=size(Ares1stack,1);

Ares1 = zeros(BackupExp.imgNumRows, BackupExp.imgNumCols, BackupExp.numFreqs);

%crop
Ares1 = Ares1(SpanYTrans,SpanXTrans,:);

for i = 1:BackupExp.numFreqs
    Ares1(:,:,i) = transpose( reshape(Ares1stack(i, :), [BackupExp.imgNumCols, BackupExp.imgNumRows] )  );
end

sizeXY = size(BinImage(Ares1(:,:,1),binSize));
bAres1 = zeros(sizeXY(1),sizeXY(2),length(Fres1)); %Bres1 = Ares1;

for i = 1:length(Fres1)
    bAres1(:,:,i) = BinImage(Ares1(:,:,i),binSize);
end

sizeX = size(bAres1,2); sizeY = size(bAres1,1); %Image dimensions
    

% Correct for severely non-unity baseline by dividing pixelwise by
%   average of all frequency points
nAres1 = zeros(size(bAres1));
NormalizationFactor1 = mean(bAres1,3);    % compute average
for i = 1:length(Fres1)
    nAres1(:,:,i) = bAres1(:,:,i) ./ NormalizationFactor1;
end

%Find global averaged spectra
spec1data = squeeze(mean(mean(nAres1,1),2));%global spectrum 51x1 array
leftbaseline=mean(spec1data(1:5,1));
rightbaseline=mean(spec1data(size(spec1data,1)-5:end,1));
globalmeanbaseline=mean([leftbaseline rightbaseline]);
baselinerange=globalmeanbaseline-min(spec1data);%difference between baseline and the minimum value in global spectrum
spec1dataZBL=spec1data-globalmeanbaseline; %ZBL=zero baseline
gbAres1temp=zeros(size(nAres1));
for i=1:size(nAres1,1)
    for j=1:size(nAres1,2)
        leftbaseline=mean(nAres1(i,j,1:5));
        rightbaseline=mean(nAres1(i,j,size(spec1data,1)-5:end));
        meanbaseline=mean([leftbaseline rightbaseline]);
        pixelrange=meanbaseline-min(nAres1(i,j,:));%difference between baseline and the minimum value in pixel spectrum
        for k=1:size(nAres1,3)
            gbAres1temp(i,j,k)=(baselinerange / pixelrange) * (nAres1(i,j,k) - meanbaseline) - 0.5*spec1dataZBL(k,1) + globalmeanbaseline;
        end
    end
end
