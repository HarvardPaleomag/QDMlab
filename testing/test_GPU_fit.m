%simple GPUfit example for testing purposes
model_id = ModelID.GAUSS_1D;
data = load("D:\data\mike\MSM16667a\P6\run_00000.mat");

%%
data_ =  data.imgStack1(:,1:100);
ig = get_initial_guess(data_, data.freqList(1:101));
%%
[initialGuess, states, chiSquares, n_iterations, time] = gpufit(data_, [], ...
            model_id, ig, ...
            [], EstimatorID.MLE, freq_);
        
function initialGuess = get_initial_guess(gpudata, freq, diamond)
arguments
    gpudata
    freq
    diamond = 'DAC'
end
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
