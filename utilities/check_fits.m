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
dataPosLeft = reshapeImg(rawDataPos.imgStack1, rawDataPos.numFreqs, rawDataPos.imgNumCols, rawDataPos.imgNumRows);
dataPosRight = reshapeImg(rawDataPos.imgStack2, rawDataPos.numFreqs, rawDataPos.imgNumCols, rawDataPos.imgNumRows);
dataNegLeft = reshapeImg(rawDataNeg.imgStack1, rawDataNeg.numFreqs, rawDataNeg.imgNumCols, rawDataNeg.imgNumRows);
dataNegRight = reshapeImg(rawDataNeg.imgStack2, rawDataNeg.numFreqs, rawDataNeg.imgNumCols, rawDataNeg.imgNumRows);
freqList = reshape(rawDataPos.freqList, [rawDataPos.numFreqs, 2]);

% prefilter data for hot pixels
[bool, dataName, ledName] = is_B111(finalFits);
data = filter_hot_pixels(finalFits.(dataName));
binning = detect_binning(finalFits);

LED = finalFits.(ledName);
ratio = min(size(LED))/max(size(LED));
% Create image
f = figure('Units', 'normalized');
set(gcf,'OuterPosition',[0.15,0.25,0.7,0.4*ratio]);

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
        x = click(1);
        y = click(2);
        titleTxt = ['X:' num2str(round(x)) '(' num2str(round(x) * binning) ')' ...
                   ' Y:' num2str(round(y)) '(' num2str(round(y) * binning) ')'];
        ax1 = subplot(1,3,1);
        
        title(titleTxt)
          
        if point ~= 0
            delete(point);
        end
        
        hold on
        px = round(x*binning);
        py = round(y*binning);
        
        point = scatter(round(x),round(y),'k');
        
        ax2 = subplot(1,3,2);
        cla()
        title(titleTxt)

        hold on
        f1 = freqList(:,1);
        idx = int16(xy2index(py/binning,px/binning, size(finalFits.negDiff)));
        plot(f1, dataPosLeft(:,px,py))
        model = model_GPU(finalFits.leftPos.p(:,idx), f1,'diamond', finalFits.kwargs.diamond)...
            / max(model_GPU(finalFits.leftPos.p(:,idx), f1, 'diamond', finalFits.kwargs.diamond))...
            * max(dataPosLeft(:,px,py));
        plot(f1, model)
        plot(f1, dataNegLeft(:,px,py))
        
        legend('+','-');
        ylabel('Intensity')
        xlabel('f (Hz)')

        ax3 = subplot(1,3,3);
        cla()
        title(titleTxt)

        hold on
        plot(freqList(:,2), dataPosRight(:,round(x*binning),round(y*binning)))
        plot(freqList(:,2), dataNegRight(:,round(x*binning),round(y*binning)))
        ylabel('Intensity')
        xlabel('f (Hz)')
    end
end