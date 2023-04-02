function [imOut] = QuadBGsub(imIn)
%[imOut] = QuadBGsub(imIn)
% Perform quadratic background subtraction on an input image.
%
% INPUTS:
% - imIn: A 2D array representing the input image. The array should have
%   dimensions (Ny, Nx), where Ny is the number of rows and Nx is the
%   number of columns.
%
% OUTPUTS:
% - imOut: A 2D array representing the input image with the quadratic
%   background subtracted out. The array has the same dimensions as the
%   input array imIn.
%
% NOTES:
% - This function fits a 2D quadratic polynomial to the input image, and
%   subtracts the fitted surface from the input image to remove the
%   background.
% - The function assumes that the background varies smoothly across the
%   input image and can be well-approximated by a quadratic polynomial.
% - The function uses the Curve Fitting Toolbox in MATLAB to perform the
%   polynomial fit.

% Convert input image to double precision
imIn = double(imIn);

% Generate x and y coordinates for the input image
x = 1:size(imIn, 2);
y = 1:size(imIn, 1);
[X, Y] = meshgrid(x, y);

% Prepare data for polynomial fitting
[xData, yData, zData] = prepareSurfaceData(X, Y, imIn);

% Fit 2D quadratic polynomial to the input image
[fitout, gof] = fit([xData, yData], zData, 'poly22');
cvals = coeffvalues(fitout);
fitFunction = cvals(1) + cvals(2) * X + cvals(3) * Y + ...
    cvals(4) * X .* X + cvals(5) * X .* Y + cvals(6) * Y .* Y;
    
% Subtract fitted surface from input image to obtain output image
imOut = (imIn - fitFunction);

end
