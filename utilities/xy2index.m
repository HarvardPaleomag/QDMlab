function index = xy2index(row, col, shape, kwargs)
%[index] = xy2index(row, col, shape; 'type')
% returns the index of a pixel in the gpudata array from given x,y of the
% pixel
%
% x: int
%   x(column) coordinate of pixel
% y: int
%   y(row) coordinate of pixel
% nRows: int
%    number of Rows in array
% type: str (gpu)
%    can be 'gpu' or 'binDataNorm'

arguments
    row
    col
    shape
    kwargs.type = 'gpu';
end

% if strcmp(kwargs.type, 'binDataNorm')
%     col_ = row;
%     row = col ;
%     col = col_;
% end

index = (row-1) * shape(2) + col;
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
