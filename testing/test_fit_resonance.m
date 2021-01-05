expData = load('D:\data\mike\NRM\run_00000.mat');
%% create data like fit_resonance
freq = expData.freqList(1:expData.numFreqs) / 1E9;   %everything is in GHz
n = 1; binSize = 4;
dataStack = expData.(sprintf('imgStack%i',n));
% data preparation
% X/Y of unbinned data
% Note: X = COL; Y = ROW -> (y,x) for (row, col) matlab convention
spanXTrans = 1:expData.imgNumCols;
spanYTrans = 1:expData.imgNumRows;

% check for 101 frequencies. File includes imgStack3
if isfield(expData, 'imgStack3')      
    % combine 1&2 or 3&4
    dataStacka = expData.(sprintf('imgStack%i',n)); 
    dataStackb = expData.(sprintf('imgStack%i',n+1));
    dataStack = [dataStacka; dataStackb];
end

data = zeros(expData.imgNumRows, expData.imgNumCols, expData.numFreqs);

% crop
data = data(spanYTrans,spanXTrans,:);

% reshape and transpose each image
for y = 1:expData.numFreqs
    data(:,:,y) = transpose(reshape(dataStack(y, :), [expData.imgNumCols, expData.imgNumRows]));
end

% binning
fprintf('<>   %i: binning data >> binSize = %i\n', n, binSize);

sizeXY = size(BinImage(data(:,:,1),binSize));
binData = zeros(sizeXY(1),sizeXY(2),length(freq));

for y = 1:length(freq)
    binData(:,:,y) = BinImage(data(:,:,y),binSize);
end

sizeX = size(binData,2); % binned image x-dimensions
sizeY = size(binData,1); % binned image y-dimensions

% Correct for severely non-unity baseline by dividing pixelwise by
% average of all frequency points

binDataNorm = zeros(size(binData));
NormalizationFactor = mean(binData,3);    % compute average

for y = 1:length(freq)
    binDataNorm(:,:,y) = binData(:,:,y) ./ NormalizationFactor;
end

% global spectra subtraction
binDataNorm = correct_global(binDataNorm, 0.5);

% first determine global guess
meanData = squeeze(mean(mean(binDataNorm,1),2));

% prepare GPUfit data
sweeplength = size(dataStack,1);
imgpts = sizeX*sizeY; % number of (x,y) pixels
gpudata = reshape(binDataNorm,[imgpts,sweeplength]); % make it into 2d matrix (row after row)
gpudata = transpose(gpudata); %transpose to make it 51 x pixels
gpudata = single(gpudata);
xValues = single(freq');
%%
disp('GLOBAL')
[fitOld, guessOld, badPixelsOld] = fit_resonance(expData, binSize, freq, n, 'type',0);
%%
disp('LOCAL WITH GLOBAL SUBSTITUTION')
[fitOldGlob, guessOldGlob, badPixelsOldGlob] = fit_resonance(expData, binSize, freq, n, 'type',1);
%%
disp('LOCAL WITH GAUSSIAN SUBSTITUTION')
[fitOldGaus, guessOldGaus, badPixelsOldGaus] = fit_resonance(expData, binSize, freq, n, 'type',1, 'gaussianFit',1);
%%
disp('LOCAL GAUSSIAN')
[fitNew, guessNew, preParamsNew] = fit_resonance(expData, binSize, freq, n, 'type',2);
%%
test(fitNew, fitOld, gpudata, binDataNorm, freq, binSize)
%%
test(fitNew, fitOldGlob, gpudata, binDataNorm, freq, binSize)
%%
test(fitNew, fitOldGaus, gpudata, binDataNorm, freq, binSize)
%%
function test(fitNew, fitOld, gpudata, binDataNorm, freq, binSize)
    close all
    x = 273;
    y = 211;
    figure
    set(gcf,'position',[150,150,1000,800])

    ax1 = subplot(3,2,1);
    old = fitOld.resonance;
%     old = filter_hot_pixels(old, 'cutOff', 20,'chi', fitOld.chiSquares, 'winSize',nan);
%     old = filter_hot_pixels(old, 'winSize',nan);
    imagesc(old,'Parent', ax1,'CDataMapping','scaled','hittest', 'off')
    title(ax1, 'OLD')

    axis equal, axis tight, axis xy
    ax1.CLim = [2.83 2.85];

    colorbar()

    ax2 = subplot(3,2,2);
    new = fitNew.resonance;
%     new = filter_hot_pixels(new, 'cutOff', 20,'chi', fitNew.chiSquares, 'winSize',nan);
%     new = filter_hot_pixels(new, 'winSize',nan);
    imagesc(new,'Parent',ax2,'CDataMapping','scaled','hittest', 'off')
    axis equal, axis tight, axis xy
    title('New')
    ax2.CLim = [2.83 2.85];
    colorbar()

    ax3 = subplot(3,2,3);
    chimax = double(max([max(fitOld.chiSquares,[],'all'), max(fitNew.chiSquares,[],'all')]));
    im = imagesc(fitOld.chiSquares,'Parent',ax3,'CDataMapping','scaled','hittest', 'off');
    axis equal, axis tight, axis xy
    title('X^2')
    ax3.CLim = [0 chimax]; 
    set(ax3,'ColorScale','log')
    colorbar()

    ax4 = subplot(3,2,4);
    imagesc(fitNew.chiSquares,'Parent',ax4,'CDataMapping','scaled','hittest', 'off')
    axis equal, axis tight, axis xy
    ax4.CLim = [0 chimax]; 
    title('X^2')
    set(ax4,'ColorScale','log')

    colorbar()
    % linkaxes([ax1 ax2 ax3 ax4])

    n = 0;
    points = [0 0 0 0 0];

    ax5 = subplot(3,2,5);
    imagesc(fitNew.resonance-fitOld.resonance,'Parent',ax5,'CDataMapping','scaled','hittest', 'off')
    mx = max(abs(fitNew.resonance-fitOld.resonance),[],'all');
    set(ax5,'ColorScale','log')
    
    colorbar()
    axis equal, axis tight, axis xy
    
    
    ax6 = subplot(3,2,6);
%     set(gcf,'WindowButtonDownFcn',@clickFitFcn)
    linkaxes([ax1 ax2 ax3 ax4 ax5]);
    for ax = [ax1 ax2 ax3 ax4 ax5]
       set(ax,'ButtonDownFcn',@clickFitFcn)
    end

    drawnow;
    function clickFitFcn(hObj, event)
        % Get click coordinate
        click = event.IntersectionPoint;
        x = round(click(1));
        y = round(click(2));

        titleTxt = ['X:' num2str(round(x)) '(' num2str(round(x) * binSize) ')' ...
                   ' Y:' num2str(round(y)) '(' num2str(round(y) * binSize) ')'];
        
        ax6 = subplot(3,2,6);
        cla()
        rows = size(binDataNorm,1);
        plot(freq, squeeze(binDataNorm(y,x,:)), 'k','DisplayName','data')
        hold on
        idx = xy2index(x, y,rows);
        plot(freq, 1+model_GPU(fitOld.p(:,idx), freq), 'r--','DisplayName','OLD fit')
        plot(freq, 1+model_GPU(fitNew.p(:,idx), freq), 'b','DisplayName','NewFit')
        plot(freq, 1+model_GPU(fitOld.g(:,idx), freq), 'r:','DisplayName','OLD IG')
        plot(freq, 1+model_GPU(fitNew.g(:,idx), freq), 'b:','DisplayName','New IG')
        ylabel('Intensity')
        xlabel('f (Hz)')
        legend('Location','bestoutside')
        for i = 1:5
            point = points(i);
            if point ~= 0
                delete(point);
            end
            subplot(3,2,i)
            hold on
            point = scatter(round(x),round(y),'xr');
            points(i) = point;
        end
        title(ax6, titleTxt)
        fprintf('%s OLD: %.5f, %.2e\t NEW: %.5f, %.2e\n', titleTxt, ...
            fitOld.resonance(y,x), fitOld.chiSquares(y,x), fitNew.resonance(y,x), fitNew.chiSquares(y,x))
    end
end
