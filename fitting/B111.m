function d = B111(finalFits, plots)
%[d] = B111(finalFits; 'resonanceFieldPlot', 'centerShiftPlot', 'fieldsPlot', 'centerShiftDistribution')
    arguments
        finalFits
        plots.resonanceFieldPlot {mustBeBoolean(plots.resonanceFieldPlot)} = false
        plots.centerShiftPlot {mustBeBoolean(plots.centerShiftPlot)} = false
        plots.fieldsPlot {mustBeBoolean(plots.fieldsPlot)} = false
        plots.centerShiftDistribution {mustBeBoolean(plots.centerShiftDistribution)} = false
    end
    
    gamma = 0.0028;
    zfs = 2.870;
    
    d.resPos1 = finalFits.leftPos.resonance;
    d.resPos2 = finalFits.rightPos.resonance;
    
    d.resNeg1 = finalFits.leftNeg.resonance;
    d.resNeg2 = finalFits.rightNeg.resonance;
    
    medianResPos = median((d.resPos2+d.resPos1)/2, 'all', 'omitnan');
    medianResNeg = median((d.resNeg2+d.resNeg1)/2, 'all', 'omitnan');
    d.med = mean([medianResPos medianResNeg]);
    
    %% shift away from the mean
    % negative fields
    d.dNeg1 = (zfs-d.resNeg1);
    d.dNeg2 = (d.resNeg2-zfs);
    % positive fields
    d.dPos1 = (zfs-d.resPos1);
    d.dPos2 = (d.resPos2-zfs);
    
    %% mean shift (i.e. mean(field) with pos/neg field applied
    d.meanPos =  (d.dPos1+d.dPos2)/2 / gamma;
    d.meanNeg =  (d.dNeg1+d.dNeg2)/2 / gamma;

    %% ferro / para component calculation
    d.ferro = (d.meanPos - d.meanNeg)/2; 
    d.para  = (d.meanPos + d.meanNeg)/2;
    
    %% center shift in pos/neg field direction
    d.deltaPos = (d.dPos2-d.dPos1)/2;
    d.deltaNeg = (d.dNeg2-d.dNeg1)/2;
    
    %% mean center shift
    d.centerShift = (d.deltaPos + d.deltaNeg)/2;
    
    %% differences between pos/ neg field
    d.shiftDelta = (d.deltaNeg - d.deltaPos)/2;
    
    %% PLOTTING
    %% center shift
    if isequal(plots.centerShiftPlot, true)
        f = figure('units','normalized','outerposition',[0.15 0.4 0.7 0.4],'NumberTitle', 'off', 'Name', 'Center Shifts');
        ax1 = subplot(1,2,1);
        QDM_figure(1000*d.deltaPos, 'ax', ax1, 'title', sprintf('Center shift (+)'), ...
            'cbTitle', '\Delta (MHz)', 'unit', 'none');
        ax2 = subplot(1,2,2);
        QDM_figure(1000*d.deltaNeg, 'ax', ax2, 'title', sprintf('Center shift (-)'), ...
            'cbTitle', '\Delta (MHz)', 'unit', 'none');
        linkaxes([ax1 ax2]);
    end
    
    %% fields
    if isequal(plots.fieldsPlot, true)
        f = figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8],'NumberTitle', 'off', 'Name', 'B111ferro, B111para');
        ax1 = subplot(2,2,1);
        QDM_figure(d.ferro, 'st', 10, 'preThreshold', 10, 'ax', ax1, 'title', 'B_{111} ferro');
        ax2 = subplot(2,2,2);
        QDM_figure(d.para-median(d.para, 'all'), 'st', 10, 'preThreshold', 10, 'ax', ax2, 'title', 'B_{111} para', 'mustBe', false);
        ax3 = subplot(2,2,3);
        QDM_figure(1000.*d.centerShift, 'st', 2, 'preThreshold', 0, 'ax', ax3, ...
            'title', 'mean(center shift)', 'unit', 'none', 'cbTitle', 'CS (MHz)', ...
            'method', 'fit', 'nOutlier', 10, 'mustBe','pos');
        ax4 = subplot(2,2,4);
        QDM_figure(d.shiftDelta, 'st', 2, 'preThreshold', 1, 'ax', ax4, 'title', 'CS^- - CS^+');
        linkaxes([ax1 ax2 ax3 ax4]);
    end
    
    %% resonance fields
    if isequal(plots.resonanceFieldPlot, true)
        f = figure('units','normalized','outerposition',[0.15 0.2 0.6 0.75],'NumberTitle', 'off', 'Name', 'Resonance Fields');
        ax1 = subplot(2,2,1);
        QDM_figure(d.dPos1-median(d.dPos1,'all'), 'ax', ax1, 'title', 'B^+_1');
        ax2 = subplot(2,2,2);
        QDM_figure(d.dPos2-median(d.dPos2,'all'), 'ax', ax2, 'title', 'B^+_2');
        ax3 = subplot(2,2,3);
        QDM_figure(d.dNeg1-median(d.dNeg1,'all'), 'ax', ax3, 'title', 'B^-_1');
        ax4 = subplot(2,2,4);
        QDM_figure(d.dNeg2-median(d.dNeg2,'all'), 'ax', ax4, 'title', 'B^-_2');
        linkaxes([ax1 ax2 ax3 ax4]);
    end
    
    if isequal(plots.centerShiftDistribution, true)
        f = figure('units','normalized','outerposition',[0.15 0.2 0.4 0.5],'NumberTitle', 'off', 'Name', 'centerShift dist.');
        a1 = subplot(1,2,1);
        a2 = subplot(1,2,2);
        
        ferro = QuadBGsub(d.ferro);
        ferro = filter_hot_pixels(ferro,'cutoff',2, 'win', nan);
        para = QuadBGsub(d.para);
        para = filter_hot_pixels(para,'cutoff',2, 'win', nan);
        CS = QuadBGsub(d.centerShift);
        CS = filter_hot_pixels(CS,'cutoff',2, 'win', nan);
        
        %% ferro
        ferroData = [reshape(ferro,[],1), reshape(CS,[],1)];
        ferroData = rmmissing(ferroData);
        
        hist3(a1, ferroData, 'Nbins',[100,100], 'CdataMode','auto','EdgeColor', 'none')
        xlabel(a1,'ferro (2\sigma filtered)')
        ylabel(a1,'centerShift (2\sigma filtered)')
        view(a1, 2)
        colormap(a1, 'turbo')
        
        %% para
        paraData = [reshape(para,[],1), reshape(CS,[],1)];
        paraData = rmmissing(paraData);
        
        hist3(a2, paraData, 'Nbins',[100,100], 'CdataMode','auto','EdgeColor', 'none')
        xlabel(a2,'para (2\sigma filtered)')
        ylabel(a2,'centerShift (2\sigma filtered)')
        colormap(a2, 'turbo')
        view(a2, 2)
    end
end