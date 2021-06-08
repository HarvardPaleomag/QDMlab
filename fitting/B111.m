function d = B111(expData, plots)
    arguments
        expData
        plots.resonanceFieldPlot {mustBeBoolean(plots.resonanceFieldPlot)} = false
        plots.centerShiftPlot {mustBeBoolean(plots.centerShiftPlot)} = false
        plots.fieldsPlot {mustBeBoolean(plots.fieldsPlot)} = false
        plots.centerShiftDistribution {mustBeBoolean(plots.centerShiftDistribution)} = false
    end
    
    gamma = 0.0028;
    zfs = 2.870;
    
    d.resPos1 = expData.leftPos.resonance;
    d.resPos2 = expData.rightPos.resonance;
    
    d.resNeg1 = expData.leftNeg.resonance;
    d.resNeg2 = expData.rightNeg.resonance;
    
    medianResPos = median((d.resPos2+d.resPos1)/2, 'all', 'omitnan');
    medianResNeg = median((d.resNeg2+d.resNeg1)/2, 'all', 'omitnan');
    d.med = mean([medianResPos medianResNeg]);
    
    %% shift away from the mean
    % negative fields
    d.dNeg1 = (zfs-d.resNeg1)/gamma;
    d.dNeg2 = (d.resNeg2-zfs)/gamma;
    % positive fields
    d.dPos1 = (zfs-d.resPos1)/gamma;
    d.dPos2 = (d.resPos2-zfs)/gamma;
    
    %% mean shift (i.e. mean(field) with pos/neg field applied
    d.meanPos =  (d.dPos1+d.dPos2)/2;
    d.meanNeg =  (d.dNeg1+d.dNeg2)/2;

    %% ferro / para component calculation
    d.ferro = (d.meanPos - d.meanNeg)/2; 
    d.para  = (d.meanPos + d.meanNeg)/2;
    
    %% center shift in pos/neg field direction
    d.deltaPos = (d.dPos2-d.dPos1)/2;
    d.deltaNeg = (d.dNeg2-d.dNeg1)/2;
    
    %% average center shift
    d.centerShift = (d.deltaPos + d.deltaNeg)/2;
    
    %% differences between pos/ neg field
    d.shiftDelta = (d.deltaNeg - d.deltaPos)/2;
    
    %% PLOTTING
    %% center shift
    if isequal(plots.centerShiftPlot, true)
        f = figure('units','normalized','outerposition',[0.15 0.4 0.7 0.4],'NumberTitle', 'off', 'Name', 'Center Shifts');
        ax1 = subplot(1,2,1);
        QDM_figure(d.deltaPos-median(d.deltaPos, 'all'), 'ax', ax1, 'title', sprintf('B^+ - median (%.2e)', median(d.deltaPos, 'all')));
        ax2 = subplot(1,2,2);
        QDM_figure(d.deltaNeg-median(d.deltaNeg, 'all'), 'ax', ax2, 'title', sprintf('B^- - median (%.2e)', median(d.deltaNeg, 'all')));
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
        QDM_figure(d.centerShift, 'st', 20, 'preThreshold', 1, 'ax', ax3, 'title', 'mean(center shift)');
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