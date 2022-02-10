function ODMR_to_B111_new(folder)

data = load_ODMR_data(folder);

end

function initialGuess = get_initial_guess(gpudata, freq, diamond)
%[initialGuess] = get_initial_guess(gpudata, freq, diamond)
initialGuess = zeros(4, size(gpudata, 2), 'single');


% amplitude
mx = nanmax(gpudata);
mn = nanmin(gpudata);
initialGuess(1, :) = -abs(((mx - mn)./mx));

% center frequency
[~, idx] = sort(gpudata);
l = 7; % lowest n values
mxidx = max(idx(1:l, :));
mnidx = min(idx(1:l, :));

if strcmp(diamond, 'N15') |  strcmp(diamond, 'DAC')
    cIdx = int16((mxidx+mnidx)/2);
else
    cIdx = int16(mean(cat(1, mxidx, mnidx)));
end

center = freq(cIdx);
initialGuess(2, :) = center;

% width
if strcmp(diamond, 'DAC')
    initialGuess(3, :) = 0.004;
else
    initialGuess(3, :) = 0.0005;
end
    
% offset
initialGuess(4, :) = mx;
end
