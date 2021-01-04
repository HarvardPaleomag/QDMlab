function g = model_gauss(x, p)
% 1: Amplitude
% 2: Center
% 3: Width
% 4: Shift
    g = p(1) * exp(-(x - p(2)).^2 / (2 * p(3)^2)) + p(4);
end