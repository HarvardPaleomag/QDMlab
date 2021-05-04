function g = model_gauss(x, p)
% 1: Amplitude
% 2: Center
% 3: Width
% 4: Shift
    msg = sprintf('calculating gaussian with amplitude = %.2e, center = %.2e, width = %.2e, shift = %.2e', p(1),P(2),P(3),P(4));
    logMsg('debug',msg,1,0);
    g = p(1) * exp(-(x - p(2)).^2 / (2 * p(3)^2)) + p(4);
end