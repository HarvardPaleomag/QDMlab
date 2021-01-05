expData = load('D:\data\mike\NRM\run_00000.mat');
%% create data like fit_resonance
freq = expData.freqList(1:expData.numFreqs) / 1E9;   %everything is in GHz
n = 1; binSize = 2;
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

fprintf('<>   %i: starting parameter estimation\n', n);

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
[fitOld, guessOld, badPixelsOld] = fit_resonance(expData, binSize, freq, n, 'type',0);
[fitOldGlob, guessOldGlob, badPixelsOldGlob] = fit_resonance(expData, binSize, freq, n, 'type',1);
[fitOldGaus, guessOldGaus, badPixelsOldGaus] = fit_resonance(expData, binSize, freq, n, 'type',1, 'gaussianFit',1);
[fitNew, guessNew, badPixelsNew] = fit_resonance(expData, 1, freq, n, 'type',2);

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
    set(gcf,'position',[150,150,1500,800])

    ax1 = subplot(3,2,1);
    title('OLD')
    old = filter_hot_pixels(fitOld.resonance, 'cutOff', 20,'chi', fitOld.chiSquares, 'winSize',nan);
%     old = filter_hot_pixels(old, 'winSize',nan);
    imagesc(old,'CDataMapping','scaled');
    axis equal, axis tight, axis xy
    ax1.CLim = [2.83 2.85];

    colorbar()

    ax2 = subplot(3,2,2);
    title('New')
    new = filter_hot_pixels(fitNew.resonance, 'cutOff', 20,'chi', fitNew.chiSquares, 'winSize',nan);
%     new = filter_hot_pixels(new, 'winSize',nan);
    imagesc(new,'Parent',ax2,'CDataMapping','scaled')
    axis equal, axis tight, axis xy
    ax2.CLim = [2.83 2.85];
    colorbar()

    ax3 = subplot(3,2,3);
    title('X2')
    chimax = double(max([max(fitOld.chiSquares,[],'all'), max(fitNew.chiSquares,[],'all')]));
    imagesc(fitOld.chiSquares,'Parent',ax3,'CDataMapping','scaled')
    axis equal, axis tight, axis xy
%     ax3.CLim = [0 chimax];
    colorbar()
    
    ax4 = subplot(3,2,4);
    title('X2')
    imagesc(fitNew.chiSquares,'Parent',ax4,'CDataMapping','scaled')
    axis equal, axis tight, axis xy
%     ax4.CLim = [0 chimax];
    colorbar()
    % linkaxes([ax1 ax2 ax3 ax4])

    n = 0;
    points = [0 0 0 0 0];

    % set(ax1,'ButtonDownFcn',@clickFitFcn)
    set(ax1,'ButtonDownFcn',@clickFitFcn)


    ax5 = subplot(3,2,5);
    imagesc(fitNew.resonance-fitOld.resonance,'Parent',ax5,'CDataMapping','scaled')
    mx = max(abs(fitNew.resonance-fitOld.resonance),[],'all');
    set(ax5,'ColorScale','log')
    colorbar()
    axis equal, axis tight, axis xy
    
    
    ax6 = subplot(3,2,6);
    set(gcf,'WindowButtonDownFcn',@clickFitFcn)
    linkaxes([ax1 ax2 ax3 ax4 ax5]);

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
        plot(freq, 1+model_GPU(fitOld.g(:,idx), freq), 'r:','DisplayName','')
        plot(freq, 1+model_GPU(fitNew.g(:,idx), freq), 'b:','DisplayName','')
        ylabel('Intensity')
        xlabel('f (Hz)')
        legend()
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

        fprintf('%s OLD: %.5f, %.2e\t NEW: %.5f, %.2e\n', titleTxt, ...
            fitOld.resonance(y,x), fitOld.chiSquares(y,x), fitNew.resonance(y,x), fitNew.chiSquares(y,x))
    end
end
