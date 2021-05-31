function [x,y] = index2xy(index, nRows, kwargs)
%[x, y] = index2xy(index, nRows; 'type')
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
    nRows
    kwargs.type = 'gpu';
end

x = fix(index / nRows)+1;
y = mod(index, nRows);

if y == 0 
    y = nRows;
    x = x-1;
end

if strcmp(kwargs.type, 'binDataNorm')
    x_ = y;
    y = x ;
    x = x_;
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
