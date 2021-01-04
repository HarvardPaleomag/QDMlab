function data = reshape_QDM_data(data)
% reshapes data to freq:row:col
    data.imgStack1 = reshape(data.imgStack1, [data.numFreqs, data.imgNumCols, data.imgNumRows]);
	data.imgStack2 = reshape(data.imgStack2, [data.numFreqs, data.imgNumCols, data.imgNumRows]);