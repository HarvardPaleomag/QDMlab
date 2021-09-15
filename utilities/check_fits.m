function check_fits(Bzdata, rawDataPos, rawDataNeg)
%check_fits(Bzdata, rawDataPos, rawDataNeg)
    close all
    disp(['<>   reshaping data: ' num2str(rawDataPos.numFreqs) 'x' num2str(rawDataPos.imgNumCols) 'x' num2str(rawDataPos.imgNumRows)])
    dataPosLeft = reshapeImg(rawDataPos.imgStack1, rawDataPos.numFreqs, rawDataPos.imgNumCols, rawDataPos.imgNumRows);
    dataPosRight = reshapeImg(rawDataPos.imgStack2, rawDataPos.numFreqs, rawDataPos.imgNumCols, rawDataPos.imgNumRows);
    dataNegLeft = reshapeImg(rawDataNeg.imgStack1, rawDataNeg.numFreqs, rawDataNeg.imgNumCols, rawDataNeg.imgNumRows);
    dataNegRight = reshapeImg(rawDataNeg.imgStack2, rawDataNeg.numFreqs, rawDataNeg.imgNumCols, rawDataNeg.imgNumRows);
    freqList = reshape(rawDataPos.freqList, [rawDataPos.numFreqs, 2]);
    size(dataPosLeft)
    % prefilter data for hot pixels
    data = filter_hot_pixels(Bzdata.B111ferro);
    binning = detect_binning(Bzdata);
    
    % Create image
    f = figure;
    set(gcf,'position',[350,350,1400,300])
    
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
        point = scatter(round(x),round(y),'k');
        
        ax2 = subplot(1,3,2);
        cla()
        title(titleTxt)

        hold on
        plot(freqList(:,1), dataPosLeft(:,round(x*binning),round(y*binning)))
        plot(freqList(:,1), dataNegLeft(:,round(x*binning),round(y*binning)))
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

function out = reshapeImg(img, nFreq, nCol, nRow)
    out = reshape(img, [nFreq, nCol, nRow]);
end