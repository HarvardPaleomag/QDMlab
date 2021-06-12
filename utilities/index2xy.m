function [row, col] = index2xy(index, shape, kwargs)
%[row, col] = index2xy(index, shape; 'type')
% returns the x,y coordinates of a pixel given the index of the gpu/raw (kw: type) data array
% 
% Parameters
% ----------
%   index: int
%       index in reshaped array
%   nRows: int
%       number of Rows in unbinned array
%   type: str [gpu]
%       can be 'gpu' or 'raw'
%
% Returns
% -------
%     x: int
%         x(column) coordinate of pixel
%     y: int
%         y(row) coordinate of pixel

% See also xy2index

arguments
    index
    shape
    kwargs.type = 'gpu';
end

if strcmp(kwargs.type, 'gpu')
    col = floor(index / shape(1))+1;
    row = mod(index, shape(1));

    if row == 0 
        row = shape(1);
        col = col-1;
    end
else
    row = floor(index / shape(2))+1;
    col = mod(index, shape(2));

    if col == 0 
        col = shape(2);
        row = row-1;
    end
end

msg = sprintf('idx: %i -> row: %i, col: %i for shape(%i, %i)', index, row, col, shape(1), shape(2));
logMsg('debug',msg,1,0);
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
