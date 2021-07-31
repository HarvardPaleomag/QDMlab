function fig = gpu_fit_checkPlot(fit, binDataNorm, freq, binSize, diamond)
%[fig] = gpu_fit_checkPlot(fit, binDataNorm, freq, binSize, diamond)
    fig = figure('units','normalized','outerposition',[0 0.5 1 0.35]);
    
    %% RESONANCE
    ax1 = subplot(1,3,1);
    
    freq = fit.freq;
    binSize = fit.binSize;
    
    fitData = fit.resonance;
    fitData(fit.states~=0) = nan;

    res = imagesc(fitData,'Parent', ax1,'CDataMapping','scaled','hittest', 'off');
    set(res,'AlphaData',~isnan(fitData))
    title(ax1, 'Resonance')
    axis equal, axis tight, axis xy
    xlabel('pixel');
    ylabel('pixel');
    colorbar()
    
    %% CHI
    ax2 = subplot(1,3,2);
    imagesc(fit.chiSquares,'Parent',ax2,'CDataMapping','scaled','hittest', 'off');
    axis equal, axis tight, axis xy
    title('X^2');
    xlabel('pixel');
    ylabel('pixel');
    set(ax2,'ColorScale','log');
    colorbar();

    n = 0;
    points = [0 0 0 0 0];
    linkaxes([ax1 ax2]);

    for ax = [ax1 ax2]
       set(ax,'ButtonDownFcn',@clickSelectPixel)
    end
    ax3 = subplot(1,3,3);
    drawnow;
    
    function clickSelectPixel(hObj, event)
    %clickSelectPixel(hObj, event)
        % Get click coordinate
        click = event.IntersectionPoint;
        x = round(click(1));
        y = round(click(2));

        titleTxt = sprintf('X: %4i (%4i) Y: %4i (%4i)', ...
            round(x),round(x)*binSize,round(y), round(y)*binSize);
        
        ax3 = subplot(1,3,3);
        cla()
        rows = size(binDataNorm,1);
        plot(freq, squeeze(binDataNorm(y,x,:)), 'k.','DisplayName','data')
        hold on
        idx = xy2index(y,x, size(binDataNorm));
        plot(freq, 1+model_GPU(fit.p(:,idx), freq, 'diamond', diamond), 'b','DisplayName','Fit')
        plot(freq, 1+model_GPU(fit.initialGuess.p(:,idx), freq, 'diamond', diamond), 'g--','DisplayName','initial guess')
        plot(freq, 1+model_GPU(fit.pg(:,idx), freq, 'diamond', diamond), 'r:','DisplayName','pre guess')
        ylabel('Intensity')
        xlabel('f (Hz)')
        legend('Location','southwest', 'NumColumns',3)
        
        for i = 1:2
            point = points(i);
            if point ~= 0
                delete(point);
            end
            subplot(1,3,i)
            hold on
            point = scatter(round(x),round(y),'xr');
            points(i) = point;
        end
        title(ax3, titleTxt)
        
        msg = sprintf('%s resonance: %.5f; X^2: %.2e', titleTxt, ...
                        fit.resonance(y,x), fit.chiSquares(y,x)');
        logMsg('info',msg,1,0);

        msg = sprintf('width: %.5f; contrast %.5f; state: %i',...
            fit.width(y,x), fit.contrastA(y,x), fit.states(y,x)');
        logMsg('info',msg,1,0);
    end
end