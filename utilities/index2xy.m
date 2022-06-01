function [row, col] = index2xy(index, shape, kwargs)
%[row, col] = index2xy(index, shape; 'type')
% returns the x,y coordinates of a pixel given the index of a flattened
% array.
%
% index: int
%   index in flattened array (either gpudata or expData)
% shape: int
%   shape of the output array (e.g. binDataNorm)
% type: str (gpu)
%    can be 'gpu' or 'raw'
%    if raw: index is asssumed to be from the raw (i.e. expData) array
%    if gpu: index is assumed to be from the gpu array
%
% Returns
% -------
%     x: int
%         x(column) coordinate of pixel
%     y: int
%         y(row) coordinate of pixel
%
% See also
% --------
%   xy2index

arguments
    index
    shape
    kwargs.type = 'gpu';
end

switch kwargs.type
    case 'raw'
        row = fix(index / shape(2))+1; % without remainder
        col = mod(index, shape(2));    % remainder

        if col == 0 
            col = shape(2);
            row = row-1;
        end
    case 'gpu'
        col = fix(index / shape(1))+1; % without remainder
        row = mod(index, shape(1));    % remainder

        if row == 0 
            row = shape(1);
            col = col-1;
        end
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
