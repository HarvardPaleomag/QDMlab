function check_ODMR(folder, laser, rawDataPos, rawDataNeg)
%check_ODMR(folder, laser, rawDataPos, rawDataNeg; 'laser', 'rawDataPos', 'rawDataNeg')
    
arguments
    folder
    laser = true
    rawDataPos = false
    rawDataNeg = false
end
    
    if isequal(rawDataPos, false)
        rawDataPos = load(fullfile(folder, 'run_00000.mat'));
    end
    if isequal(rawDataNeg, false)
        rawDataNeg = load(fullfile(folder, 'run_00001.mat'));
    end
    
    close all
    
    dataPosLeft = QDMreshape(rawDataPos.imgStack1, rawDataPos.imgNumRows, rawDataPos.imgNumCols);
    dataPosRight = QDMreshape(rawDataPos.imgStack2, rawDataPos.imgNumRows, rawDataPos.imgNumCols);
    dataNegLeft = QDMreshape(rawDataNeg.imgStack1, rawDataPos.imgNumRows, rawDataNeg.imgNumCols);
    dataNegRight = QDMreshape(rawDataNeg.imgStack2, rawDataPos.imgNumRows, rawDataNeg.imgNumCols);
    
    mask = any(~isnan(dataPosLeft),3) & any(~isnan(dataPosRight),3) & any(~isnan(dataNegLeft),3) & any(~isnan(dataNegRight),3);
    
    meanPosLeft = squeeze(mean(crop_data(dataPosLeft, mask), [1,2], 'omitnan'));
    meanPosRight = squeeze(mean(crop_data(dataPosRight, mask), [1,2], 'omitnan'));
    meanNegLeft = squeeze(mean(crop_data(dataNegLeft, mask), [1,2], 'omitnan'));
    meanNegRight = squeeze(mean(crop_data(dataNegRight, mask), [1,2], 'omitnan'));

    freqList = reshape(rawDataPos.freqList, [rawDataPos.numFreqs, 2]);
%     freqList = [rawDataPos.freqList;rawDataPos.freqList]';

    laser = load(fullfile(folder, 'laser.csv'));        
    led = load(fullfile(folder, 'led.csv'));
    
    % Create image
    f = figure;
    set(gcf,'position',[350,350,1400,900]);
    
    % plot QDM data
    ax1 = subplot(2,2,1);
    LED = imagesc(led,'Parent',ax1,'CDataMapping','scaled');
    colormap(ax1, 'bone');
    axis equal, axis tight, axis xy
    
    ax2 = subplot(2,2,2);
    LASER = imagesc(laser,'Parent',ax2,'CDataMapping','scaled');
    colormap(ax2, 'gray');
    linkaxis([ax1 ax2]);
    
    axis equal, axis tight, axis xy
    
    n = 0;
    points = 0;
    set(LED,'ButtonDownFcn',@buttonSelectPixel)
    set(LASER,'ButtonDownFcn',@buttonSelectPixel)
        
    points = [0, 0];
    function buttonSelectPixel(hObj, event)
    %buttonSelectPixel(hObj, event)
        % Get click coordinate
        click = event.IntersectionPoint;
        x = click(1);
        y = click(2);
        titleTxt = ['X:' num2str(round(x)) ')' ...
                   ' Y:' num2str(round(y)) ')'];
        
        for i = 1:2
            point = points(i);
            if point ~= 0
                delete(point);
            end
            subplot(2,2,i)
            hold on
            point = scatter(round(x),round(y),'xr');
            points(i) = point;
        end
        
        ax3 = subplot(2,2,3);
        cla()
        title(titleTxt)

        hold on
        plot(ax3, freqList(:,1), meanPosLeft, '.--')
        plot(ax3, freqList(:,1), meanNegLeft, '.--')
        plot(ax3, freqList(:,1), squeeze(dataPosLeft(round(y), round(x), :)))
        plot(ax3, freqList(:,1), squeeze(dataNegLeft(round(y), round(x), :)))
        legend('mean(+)', 'mean(-)', '+', '-');
        ylabel('Intensity')
        xlabel('f (Hz)')

        ax4 = subplot(2,2,4);
        cla()
        title(titleTxt)

        hold on
        plot(ax4, freqList(:,1), meanPosRight, '.--')
        plot(ax4, freqList(:,1), meanNegRight, '.--')
        plot(ax4, freqList(:,2), squeeze(dataPosRight(round(y), round(x),:)));
        plot(ax4, freqList(:,2), squeeze(dataNegRight(round(y), round(x),:)))
        legend('mean(+)', 'mean(-)', '+', '-');

        ylabel('Intensity')
        xlabel('f (Hz)')
    end
end

function croppedData = crop_data(data, mask)
%[croppedData] = crop_data(data, mask)
    [x, y, w, h] = get_mask_extent(mask);

    % cut around data
    croppedData = data(y:y+h, x:x+w,:);
end