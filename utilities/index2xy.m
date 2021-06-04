function [row, col] = index2xy(index, shape, kwargs)
%[row, col] = index2xy(index, shape; 'type')
% returns the x,y coordinates of a pixel given the index of the gpu array
%
% index: int
%   index in reshaped array
% nRows: int
%   number of Rows in unbinned array
% Returns
% -------
%     x: int
%         x(column) coordinate of pixel
%     y: int
%         y(row) coordinate of pixel
% type: str (gpu)
%    can be 'gpu' or 'binDataNorm'
% See also xy2index

arguments
    index
    shape
    kwargs.type = 'gpu';
end

row = fix(index / shape(1))+1;
col = mod(index, shape(1));

if col == 0 
    col = shape(2);
    row = row-1;
end

if strcmp(kwargs.type, 'binDataNorm')
    x_ = col;
    col = row ;
    row = x_;
end
end


%% test
% [r,c] = size(binDataNorm);
% n = 2;
% a = gpudata(:,n);
% 
% [x,y] = index2xy(n,r);
% a_ = squeeze(binDataNorm(y,x,:));
% [x,y] = index2xy(200, r)
% idx = xy2index(x,y,r)
