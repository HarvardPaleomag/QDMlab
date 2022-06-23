function index = xy2index(row, col, shape, type)
%[index] = xy2index(row, col, shape; 'type')
% returns the index of a pixel in the data (gpu/raw) array from given x,y of the
% pixel
%
% x: int
%   x(column) coordinate of pixel
% y: int
%   y(row) coordinate of pixel
% shape: (int int)
%    number of Rows, Cols in array binDataNorm array
% type: str (gpu)
%    can be 'gpu' or 'raw'. 
%    if gpu: it will return the index used in the gpu array
%    if raw: it will return the index of the raw data array

arguments
    row
    col
    shape
    type = 'gpu';
end

switch type
    case 'raw'
        index = (row-1) * shape(2) + col;
    case 'gpu'
        index = (col-1) * shape(1) + row;
end


msg = sprintf('row: %i, col: %i -> idx: %i for shape(%i, %i)', row, col, index, shape(1), shape(2));
logMsg('debug',msg,1,0);
end

%% test
% [r,c] = size(binDataNorm);
% x = 14; y =22;
% a_ = squeeze(binDataNorm(x,y,:));
% n = xy2index(x,y,r);
% a = gpudata(:,n);
% [x,y] = index2xy(200, r)
% idx = xy2index(x,y,r)
