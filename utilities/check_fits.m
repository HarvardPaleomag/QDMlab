function check_fits(finalFits, rawDataPos, rawDataNeg)
%check_fits('finalFits', 'rawDataPos', 'rawDataNeg')
% 
% Parameters
% ----------
%   finalFits:
%   rawDataPos:
%   rawDataNeg:
%   finalFits: ('none')
%   rawDataPos: ('none')
%   rawDataNeg: ('none')
% 
% Returns
% ----------

arguments
    finalFits = 'none'
    rawDataPos = 'none'
    rawDataNeg = 'none'
end

%% load all the data
finalFits = automatic_input_ui__(finalFits, 'type', 'file', ...
    "filter",'final_fits*.mat', 'title', 'select the final_fits*.mat file.', ...
    'single', 1);
rawDataPos = automatic_input_ui__(rawDataPos, 'type', 'file', ...
    "filter",'run_0000*.mat', 'title', 'select the run_00000.mat file.', ...
    'single', 1);
rawDataNeg = automatic_input_ui__(rawDataNeg, 'type', 'file', ...
    "filter",'run_0000*.mat', 'title', 'select the run_00001.mat file.', ...
    'single', 1);

if ~isstruct(finalFits); finalFits = load(finalFits); end
if ~isstruct(rawDataPos); rawDataPos = load(rawDataPos); end
if ~isstruct(rawDataNeg); rawDataNeg = load(rawDataNeg); end

%% figure
close all
disp(['<>   reshaping data: ' num2str(rawDataPos.numFreqs) 'x' num2str(rawDataPos.imgNumCols) 'x' num2str(rawDataPos.imgNumRows)])

binSize = detect_binning(finalFits);
dataPosLeft = prepare_raw_data(rawDataPos, binSize, 1);
dataPosRight = prepare_raw_data(rawDataPos, binSize, 2);
dataNegLeft = prepare_raw_data(rawDataNeg, binSize, 1);
dataNegRight = prepare_raw_data(rawDataNeg, binSize, 2);

% dataPosLeft = reshapeImg(rawDataPos.imgStack1, rawDataPos.numFreqs, rawDataPos.imgNumCols, rawDataPos.imgNumRows);
% dataPosRight = reshapeImg(rawDataPos.imgStack2, rawDataPos.numFreqs, rawDataPos.imgNumCols, rawDataPos.imgNumRows);
% dataNegLeft = reshapeImg(rawDataNeg.imgStack1, rawDataNeg.numFreqs, rawDataNeg.imgNumCols, rawDataNeg.imgNumRows);
% dataNegRight = reshapeImg(rawDataNeg.imgStack2, rawDataNeg.numFreqs, rawDataNeg.imgNumCols, rawDataNeg.imgNumRows);
freqList = reshape(rawDataPos.freqList, [rawDataPos.numFreqs, 2]);

% prefilter data for hot pixels
[~, dataName, ledName] = is_B111(finalFits);
data = filter_hot_pixels(finalFits.(dataName));
binning = detect_binning(finalFits);

LED = finalFits.(ledName);
ratio = min(size(LED))/max(size(LED));
% Create image
f = figure('Units', 'normalized');
set(gcf,'OuterPosition',[0.15,0.25,0.7,0.5*ratio]);
movegui(f,'center')

% plot QDM data
ax1 = subplot(1,3,1);
qdm = imagesc(data,'Parent',ax1,'CDataMapping','scaled');
axis equal, axis tight, axis xy

n = 0;
point = 0;
set(qdm,'ButtonDownFcn',@buttonSelectPixel)

ax2 = subplot(1,3,2);
ax3 = subplot(1,3,3);

function buttonSelectPixel(hObj, event)
    %buttonSelectPixel(hObj, event)
        % Get click coordinate
        click = event.IntersectionPoint;
        x = round(click(1));
        y = round(click(2));
        titleTxt = ['X:' num2str(round(x)) '(' num2str(round(x) * binning) ')' ...
                   ' Y:' num2str(round(y)) '(' num2str(round(y) * binning) ')'];
        ax1 = subplot(1,3,1);
        
        title(titleTxt)
          
        if point ~= 0
            delete(point);
        end
        
        hold on
        
        point = scatter(round(x),round(y),'k');
        
        ax2 = subplot(1,3,2);
        cla()
        title(titleTxt)

        hold on
        f1 = finalFits.leftPos.freq;
        f2 = finalFits.rightPos.freq;

        idx = int16(xy2index(y, x, size(dataPosRight), 'gpu'));

        dNL = squeeze(dataPosLeft( y, x, :));
        dPL = squeeze(dataNegLeft( y, x, :));

        mPL = model_GPU(finalFits.leftPos.parameters(:, y, x), f1,'diamond', finalFits.kwargs.diamond);
        mPL = normalize(mPL, 'range', [min(dPL), max(dPL)]);
        mNL = 1 + model_GPU(finalFits.leftNeg.parameters(:, y, x), f1,'diamond', finalFits.kwargs.diamond);
        mNL = normalize(mNL, 'range', [min(dNL), max(dNL)]);

        dNR = squeeze(dataPosRight(y, x, :));
        dPR = squeeze(dataNegRight(y, x, :));

        mPR = 1 + model_GPU(finalFits.rightPos.parameters(:, y, x), f2,'diamond', finalFits.kwargs.diamond);
        mPR = normalize(mPR, 'range', [min(dPR), max(dPR)]);

        mNR = 1 + model_GPU(finalFits.rightNeg.parameters(:, y, x), f2,'diamond', finalFits.kwargs.diamond);
        mNR = normalize(mNR, 'range', [min(dNR), max(dNR)]);

        plot(f1, dNL / max(dNL), 'r.-', 'DisplayName','-')
        plot(f1, mNL / max(mNL), 'r--', 'DisplayName', '- fit')

        plot(f1, dPL / max(dPL), 'b.-', 'DisplayName','+')
        plot(f1, mPL / max(mPL), 'b--', 'DisplayName', '+ fit')
        
        legend();
        ylabel('Intensity')
        xlabel('f (Hz)')

        ax3 = subplot(1,3,3);
        cla()
        hold on
        title(titleTxt)

        plot(f2, dNR / max(dNR), 'r.-', 'DisplayName','-')
        plot(f2, mNR / max(mNR), 'r--')

        plot(f2, dPR / max(dPR), 'b.-', 'DisplayName','+')
        plot(f2, mPR / max(mPR), 'b--')

        ylabel('Intensity')
        xlabel('f (Hz)')
    end
end