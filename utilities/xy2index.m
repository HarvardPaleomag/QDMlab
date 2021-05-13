function index = xy2index(x,y, nRows, kwargs)
%[index] = xy2index(x, y, nRows; 'type')
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
    x
    y
    nRows
    kwargs.type = 'gpu';
end

if strcmp(kwargs.type, 'binDataNorm')
    x_ = y;
    y = x ;
    x = x_;
end
index = (x-1)*nRows + y;

end

%% test
% [r,c] = size(binDataNorm);
% x = 14; y =22;
% a_ = squeeze(binDataNorm(x,y,:));
% n = xy2index(x,y,r);
% a = gpudata(:,n);
% [x,y] = index2xy(200, r)
% idx = xy2index(x,y,r)
