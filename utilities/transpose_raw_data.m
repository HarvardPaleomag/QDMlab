function Ares1 = transpose_raw_data(expData, nRes)
% transposes the freq x n data from the QDM into a n x m x freq large 3D image

SpanXTrans = 1:expData.imgNumCols;
SpanYTrans = 1:expData.imgNumRows;

%Data loading
freqs = expData.freqList(1:expData.numFreqs)/1E9;   %everything is in GHz

dataStack = expData.(sprintf('imgStack%i',nRes));

Ares1 = zeros(expData.imgNumRows, expData.imgNumCols, expData.numFreqs);

%crop
Ares1 = Ares1(SpanYTrans,SpanXTrans,:);

for i = 1:expData.numFreqs
    Ares1(:,:,i) = transpose( reshape(dataStack(i, :), [expData.imgNumCols, expData.imgNumRows] )  );
end