function globalFraction = globalFraction_estimator(binDataNorm)
    
    %% remove non fitting pixels
    % transform data back into picel by pixel col by col data
    nFreq = size(binDataNorm, 3);
    imgStack = QDMreshape_reverse(binDataNorm, nFreq);
    stdD = nanstd(imgStack);
    pixErr = stdD < 1e-3; % find idx of pixels with std close to zero -> all freqs have same value
       
    %% determine indices
    [~, idxMin] = nanmin(imgStack);
    idxMin(pixErr) = nan;

    [~, iMin] = nanmin(idxMin); % left most minimum
    [~, iMax] = nanmax(idxMin); % right most minimum
    
    %%
    globalMean = squeeze(mean(binDataNorm, [1,2]));
    nCol = size(binDataNorm, 2);
    
    %% figure
    % initialize UI figure
    f = figure('Position',[200 200 800 275],'NumberTitle', 'off', 'Name', 'global correction: 0.00');
    
    % define axes objects
    ax1 = axes('Parent',f,'position',[0.1 0.3  0.25 0.54]);
    ax2 = axes('Parent',f,'position',[0.4 0.3  0.25 0.54]);
    ax3 = axes('Parent',f,'position',[0.7 0.3  0.25 0.54]);

    [x1,y1] = index2xy(iMin, nCol, 'type', 'binDataNorm');
    plot(ax1, squeeze(binDataNorm(y1,x1,:)), 'k','lineWidth',1)
    hold(ax1, 'on')
    plot(ax1, globalMean, 'b:')
    p1 = plot(ax1, squeeze(binDataNorm(y1,x1,:)),'lineWidth',1);
    title(ax1, 'min(min) pixel')

    [x2,y2] = index2xy(randi(size(idxMin,2)), nCol, 'type', 'binDataNorm');
    plot(ax2, squeeze(binDataNorm(y2,x2,:)), 'k','lineWidth',1);
    hold(ax2, 'on');
    plot(ax2, globalMean, 'b:');
    p2 = plot(ax2, squeeze(binDataNorm(y2,x2,:)),'lineWidth',1);
    title(ax2, 'random pixel');

    [x3,y3] = index2xy(iMax, nCol, 'type', 'binDataNorm');
    plot(ax3, squeeze(binDataNorm(y3,x3,:)), 'k','lineWidth',1);
    hold(ax3, 'on');
    plot(ax3, globalMean, 'b:');
    p3 = plot(ax3, squeeze(binDataNorm(y3,x3,:)),'lineWidth',1);
    title(ax3, 'max(min) pixel');

    sld = uicontrol(f, 'Position',[100 50 600 3],'Style','slider','Min',0,'Max',1,'SliderStep',[0.01 0.10]);
    addlistener(sld, 'Value', 'PostSet',@(src,event) updatePlot(src, event, binDataNorm, globalMean, x1,x2,x3,y1,y2,y3, p1,p2,p3,f));
    
    sld.Value = 0;

    % Create ValueChangedFcn callback
    function updatePlot(src, event, binDataNorm, globalMean, x1,x2,x3,y1,y2,y3, p1,p2,p3, f)
        glob = sld.Value;
        dmin = correct_global(binDataNorm(y1,x1,:), glob, 'mean', globalMean);
        drand = correct_global(binDataNorm(y2,x2,:), glob, 'mean', globalMean);
        dmax = correct_global(binDataNorm(y3,x3,:), glob, 'mean', globalMean);

        set(p1,'ydata',squeeze(dmin));
        set(p2, 'ydata', squeeze(drand));
        set(p3, 'ydata', squeeze(dmax));
        set(f, 'Name', sprintf('global correction: %.3f', glob));
        drawnow;
    end
end