function fig = gpu_fit_checkPlot(fit, binDataNorm, freq, binSize)
    fig = figure;
    set(gcf,'position',[250,500,1500,250])

    ax1 = subplot(1,3,1);
    
    freq = fit.freq;
    binSize = fit.binSize;
    
    fitData = fit.resonance;
%     fitData = filter_hot_pixels(fitData, 'cutOff', 20,'chi', fitOld.chiSquares, 'winSize',nan);
%     fitData = filter_hot_pixels(fitData, 'winSize',nan);
    imagesc(fitData,'Parent', ax1,'CDataMapping','scaled','hittest', 'off')
    title(ax1, 'Resonance')
    axis equal, axis tight, axis xy
%     ax1.CLim = [2.83 2.85];

    colorbar()

    ax2 = subplot(1,3,2);
    imagesc(fit.chiSquares,'Parent',ax2,'CDataMapping','scaled','hittest', 'off');
    axis equal, axis tight, axis xy
    title('X^2')
    set(ax2,'ColorScale','log')
    colorbar()

    n = 0;
    points = [0 0 0 0 0];
    linkaxes([ax1 ax2]);

    for ax = [ax1 ax2]
       set(ax,'ButtonDownFcn',@clickFitFcn)
    end
    ax3 = subplot(1,3,3);
    drawnow;
    function clickFitFcn(hObj, event)
        % Get click coordinate
        click = event.IntersectionPoint;
        x = round(click(1));
        y = round(click(2));

        titleTxt = sprintf('X: %4i (%4i) Y: %4i (%4i)', ...
            round(x),round(x)*binSize,round(y), round(y)*binSize);
        
        ax3 = subplot(1,3,3);
        cla()
        rows = size(binDataNorm,1);
        plot(freq, squeeze(binDataNorm(y,x,:)), 'k.-','DisplayName','data')
        hold on
        idx = xy2index(x, y,rows);
        plot(freq, 1+model_GPU(fit.p(:,idx), freq), 'b','DisplayName','Fit')
        plot(freq, 1+model_GPU(fit.g(:,idx), freq), 'g--','DisplayName','initial guess')
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
        fprintf('<>            %s resonance: %.5f X^2: %.2e\n', titleTxt, ...
            fit.resonance(y,x), fit.chiSquares(y,x))
    end
end