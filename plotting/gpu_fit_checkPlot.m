<<<<<<< Updated upstream
function fig = gpu_fit_checkPlot(fit, binDataNorm, freq, binSize, diamond)
%[fig] = gpu_fit_checkPlot(fit, binDataNorm, freq, binSize, diamond)
    fig = figure('units','normalized','outerposition',[0 0.5 1 0.35]);
    
    % make str from list
    gfs = {size(GFs,2)};
    for i = 1:size(GFs,2)
        GF = GFs(i);
        gfs{i} = mat2str(GF);
    end

    fig = figure('units','normalized','outerposition',[0.05 0.5 0.9 0.3]);
    ax1 = subplot(1,3,1);
    ax2 = subplot(1,3,2);
    ax3 = subplot(1,3,3);

    c = uicontrol(fig, ...
        'Units', 'normalized', 'Position',[0.95 0.9 0.05 0.1],...
        'Style','popupmenu');
    c.String = gfs;
    c.Callback = @changeSelection;
    plot_maps(c.Value);
    points = [0 0 0 0 0];
    
    function changeSelection(src, event)
        idxGF = src.Value;
        plot_maps(idxGF)
        plot_spectra(x,y)
    end
    
    function clickSelectPixel(hObj, event)
    %clickSelectPixel(hObj, event)
        % Get click coordinate
        click = event.IntersectionPoint;

        x = round(click(1));
        y = round(click(2));
        plot_spectra(x,y);
    end

    function plot_maps(idxGF)
        d = DATA(:,:,:,idxGF);
        %% RESONANCE
        freq = fit.freq(1,:);
        binSize = fit.binSize;

        fitData = fit.resonance(:,:,idxGF);
        fitData(fit.states(:,:,idxGF)~=0) = nan;

        res = imagesc(fitData, 'Parent', ax1,'CDataMapping','scaled','hittest', 'off');
        set(res,'AlphaData',~isnan(fitData))
        title(ax1, 'Resonance')

        %% CHI
        imagesc(fit.chiSquares(:,:,idxGF),'Parent',ax2,'CDataMapping','scaled','hittest', 'off');
        title('X^2');

        for ax = [ax1 ax2]
            axis(ax, 'tight');
            axis(ax, 'equal')
            axis(ax, 'xy')
            xlabel(ax, 'pixel');
            ylabel(ax, 'pixel');
            colorbar(ax)
        end

        n = 0;
        points = [0 0 0 0 0];
        linkaxes([ax1 ax2]);

        for ax = [ax1 ax2]
           set(ax,'ButtonDownFcn',@clickSelectPixel)
        end
        drawnow;
    end
    
    function plot_spectra(x,y)
        % get data                
        idxGF = c.Value;
        d = DATA(:,:,:,idxGF);
        % determine gpudata index
        idx = xy2index(y,x, size(d), 'type', 'gpu');
        
        % clean subplot
        ax3 = subplot(1,3,3);
        cla(ax3)
        
        p = reshape(fit.p, 6, [], size(GFs,2));
        pg = reshape(fit.pg, 6, [], size(GFs,2));
        initialP = reshape(fit.initialGuess.p, 6, [], size(GFs,2));

        titleTxt = sprintf('%4i (%4i,%4i)', ...
            idx, round(x),round(y));
        
 
        plot(freq, squeeze(d(:,y,x)), 'k.','DisplayName','data', 'LineWidth', 1.5)
        hold on
        plot(freq, 1+model_GPU(p(:,idx, idxGF), freq, 'diamond', diamond), 'b','DisplayName','Fit', 'LineWidth', 1.5)
        plot(freq, 1+model_GPU(initialP(:,idx, idxGF), freq, 'diamond', diamond), 'g--','DisplayName','initial guess', 'LineWidth', 1.5)
        
        if ~isequal(fit.pg, 'none')
            plot(freq, 1+model_GPU(pg(:,idx, idxGF), freq, 'diamond', diamond), 'r:','DisplayName','pre guess', 'LineWidth', 1.5)
        end
        
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
                        fit.resonance(y,x, idxGF), fit.chiSquares(y,x, idxGF)');
        logMsg('info',msg,1,0);

        msg = sprintf('width: %.5f; contrast %.5f; state: %i',...
            fit.width(y,x), fit.contrastA(y,x, idxGF), fit.states(y,x, idxGF)');
        logMsg('info',msg,1,0);
    end
=======
function fig = gpu_fit_checkPlot(fits, data, diamond)
    %[fig] = gpu_fit_checkPlot(fit, binDataNorm, freq, binSize, diamond)
    DATA = reshape(data, [], size(fits.resonance,1), size(fits.resonance,2), size(fits.resonance,3));
    persistent x
    persistent y
    persistent ax1
    persistent ax2
    persistent ax3
    persistent ax4

    GFs = fits.GFs;
    freq = fits.freq;
    binSize = fits.binSize;
    
    % make str from list
    gfs = {size(GFs,2)};
    for i = 1:size(GFs,2)
        GF = GFs(i);
        gfs{i} = mat2str(GF);
    end
    
    %%
    compIdx = ones(size(fits.chiSquares, [1,2]));
    compGF = GFs(1)*ones(size(fits.chiSquares, [1,2]));
    for idx = 2:size(GFs,2)
        betterChi = fits.chiSquares(:,:,idx) < fits.chiSquares(:,:,idx-1);
        compIdx(betterChi) = idx;
        compGF(betterChi) = GFs(idx);
    end
    
    %%
    fig = figure('units','normalized','outerposition',[0.05 0.5 0.8 0.8]);
    ax1 = subplot(2,2,1);
    ax2 = subplot(2,2,2);
    ax3 = subplot(2,2,3);
    ax4 = subplot(2,2,4);
    
    %%
    [fitData, chiData] = get_plotData(1);
    
    res = imagesc(fitData, 'Parent', ax1,'CDataMapping','scaled','hittest', 'off');
    set(res,'AlphaData',~isnan(fitData))
    title(ax1, 'Resonance')

    %% CHI
    chi = imagesc(chiData,'Parent',ax2,'CDataMapping','scaled','hittest', 'off');
    title('X^2');

    imagesc(ax3, compGF);
    colormap(ax3, jet(size(GFs,2)));
    title(ax3,'best GF / lowest \chi^2')
    
    for a = [ax1 ax2 ax3]
        axis(a, 'tight');
        axis(a, 'equal')
        axis(a, 'xy')
        xlabel(a, 'pixel');
        ylabel(a, 'pixel');
        colorbar(a)
        set(a,'ButtonDownFcn',@clickSelectPixel);
    end
    
    set(ax4, 'Position', [ax4.Position(1:2), ax1.Position(3:4)]);
    n = 0;
    points = [0 0 0 0 0];
    linkaxes([ax1 ax2 ax3]);
    
    %%
    compIdx = ones(size(fits.chiSquares, [1,2]));
    compGF = GFs(1)*ones(size(fits.chiSquares, [1,2]));
    
    for idxGF = 2:size(GFs,2)
        betterChi = fits.chiSquares(:,:,idxGF) < fits.chiSquares(:,:,idxGF-1);
        compIdx(betterChi) = idxGF;
        compGF(betterChi) = GFs(idxGF);
    end

    c = uicontrol(fig, ...
        'Units', 'normalized', 'Position',[0.95 0.9 0.05 0.1],...
        'Style','popupmenu');
    c.String = gfs;
    c.Callback = @changeSelection;
    plot_maps(c.Value);
    points = [0 0 0 0 0 0];
    
    function changeSelection(src, event)
        idxGF = src.Value;
        plot_maps(idxGF);
        mark_selection();
        plot_spectra(x,y)
    end
    
    function clickSelectPixel(hObj, event)
    %clickSelectPixel(hObj, event)
        % Get click coordinate
        click = event.IntersectionPoint;
        x = round(click(1));
        y = round(click(2));
        mark_selection();
        plot_spectra(x,y);
        log_data();
    end
    
    function [fitData, chiData] = get_plotData(idxGF)
        d = DATA(:,:,:,idxGF);
        %% RESONANCE
        freq = fits.freq(1,:);
        binSize = fits.binSize;

        fitData = fits.resonance(:,:,idxGF);
        fitData(fits.states(:,:,idxGF)~=0) = nan;
        
        chiData = fits.chiSquares(:,:,idxGF);
    end

    function plot_maps(idxGF)
        [fitData, chiData] = get_plotData(idxGF);
        set(res,'CData', fitData);
        set(chi,'CData', chiData);
        drawnow;
    end
    
    function mark_selection()

        for axIdx = 1:3
            ax = subplot(2,2,axIdx);
            hold on
            
            % remoe old points
            point = points(axIdx);
            if point ~= 0
                delete(point);
            end
            
            % add new points
            point = scatter(ax, round(x),round(y),'Xr', 'lineWidth', 1.2);
            points(axIdx) = point;
        end
    end
    
    function log_data()
        % determine gpudata index
        d = DATA(:,:,:,idxGF);

        idx = xy2index(y,x, size(d), 'type', 'gpu');
        
        titleTxt = title_text();
        msg = sprintf('%s resonance: %.5f; X^2: %.2e', titleTxt, ...
                        fits.resonance(y,x, idxGF), fits.chiSquares(y,x, idxGF)');
        logMsg('info',msg,1,0);

        msg = sprintf('width: %.5f; contrast %.5f; state: %i',...
            fits.width(y,x), fits.contrastA(y,x, idxGF), fits.states(y,x, idxGF)');
        logMsg('info',msg,1,0);
    end

    function plot_spectra(x,y)
        % get data                
        idxGF = c.Value;
        d = DATA(:,:,:,idxGF);
        % determine gpudata index
        idx = xy2index(y,x, size(d), 'type', 'gpu');
        
        % clean subplot
        ax4 = subplot(2,2,4);
        cla(ax4)

        p = reshape(fits.p, 6, [], size(GFs,2));
        pg = reshape(fits.pg, 6, [], size(GFs,2));
        initialP = reshape(fits.initialGuess.p, 6, [], size(GFs,2));        
 
        plot(freq, squeeze(d(:,y,x)), 'k.','DisplayName','data', 'LineWidth', 1.5)
        hold(ax4, 'on')
        plot(freq, 1+model_GPU(p(:,idx, idxGF), freq, 'diamond', diamond), 'b','DisplayName','Fit', 'LineWidth', 1.5)
        plot(freq, 1+model_GPU(initialP(:,idx, idxGF), freq, 'diamond', diamond), 'g--','DisplayName','initial guess', 'LineWidth', 1.5)
        
        if ~isequal(fits.pg, 'none')
            plot(freq, 1+model_GPU(pg(:,idx, idxGF), freq, 'diamond', diamond), 'r:','DisplayName','pre guess', 'LineWidth', 1.5)
        end
        
        ylabel('Intensity')
        xlabel('f (Hz)')
        legend('Location','southwest', 'NumColumns',3)
        
        %% set title
        titleTxt = title_text();
        title(ax4, titleTxt)
    end
    function  titleTxt = title_text()
        titleTxt = sprintf('%4i (%4i,%4i)', ...
            idx, round(x),round(y));
    end
>>>>>>> Stashed changes
end