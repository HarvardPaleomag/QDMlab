function fig = gpu_fit_checkPlot(fit, data, freq, binSize, diamond, GFs)
    %[fig] = gpu_fit_checkPlot(fit, binDataNorm, freq, binSize, diamond)
    DATA = reshape(data, [], size(fit.resonance,1), size(fit.resonance,2), size(fit.resonance,3));
    persistent x
    persistent y
    
    % make str from list
    gfs = {size(GFs,2)};
    for i = 1:size(GFs,2)
        GF = GFs(i);
        gfs{i} = mat2str(GF);
    end

    fig = figure('units','normalized','outerposition',[0.05 0.5 0.9 0.33]);
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
        set(ax1, 'CLim', [min(fitData,[],'all'), max(fitData,[],'all')]);
        set(res,'AlphaData',~isnan(fitData))
        title(ax1, 'Resonance')

        %% CHI
        chiData = fit.chiSquares(:,:,idxGF);
        chi = imagesc(chiData,'Parent',ax2,'CDataMapping','scaled','hittest', 'off');
        set(ax2, 'CLim', [min(chiData,[],'all'), max(chiData,[],'all')]);
        title(ax2, 'X^2');

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

        titleTxt = sprintf('x:%4i, y:%4i (%4i) GF =%.2f ', ...
            round(x),round(y), idx, GFs(idxGF));
        
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

end