function window = return_bin_data(data, row, col, binSize)
% returns the data used for the bin binnedData(row,col,:)

% calculate the indices
rIdx = row*binSize-1:row*binSize-1+binSize-1;
cIdx = col*binSize-1:col*binSize-1+binSize-1;

% crop the data
window = data(rIdx, cIdx, :);