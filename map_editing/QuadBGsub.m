function [ imOut ] = QuadBGsub( imIn )
%Does quadratic backgound subtraction on image imIn and returns the result
msg = sprintf('subtracting quadratic background');
logMsg('debug',msg,1,0);

x = 1:size(imIn,2); y = 1:size(imIn,1);
[X,Y] = meshgrid(x,y);
[xData, yData, zData] = prepareSurfaceData( X, Y, imIn );

[fitout gof] = fit([xData, yData],zData, 'poly22');
cvals = coeffvalues(fitout);
fitFunction = cvals(1) + cvals(2)*X + cvals(3)*Y +...
     cvals(4)*X.*X + cvals(5)*X.*Y + cvals(6)*Y.*Y;
imOut = (imIn-fitFunction);

end

