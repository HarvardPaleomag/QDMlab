function filteredMaps = threshold_QDM_map(kwargs, filterArgs)
%[filteredMaps] = threshold_QDM_map('nFiles', 'checkPlot', 'save', 'removePixelAlerts', 'cutOff', 'includeHotPixel', 'chi', 'winSize', 'threshold')
% a threshold and saves it with '_thresh' tag.
%
% Parameters
% ----------
%     nFiles: cell ['none']
%     cutOff: int ['none']
%         how many standard deviations have to be exceeded for the pixel to
%         be filtered.
%     includeHotPixel: bool [false]
%         if true: the mean will be calculated including the hot pixel
%         if false: the mean is calculated after setting the pixel to nan
%     chi: array [false]
%         if chi is provided the data will be filtered according to the chi
%         values.
%     save: bool [true]
%         Specifies if the data should be saved as a new file.
%     winSize: int [3]
%         specifies the number of pixels to the left AND right to be used for
%         averaging
%         if nan: values are replaced by nan
%     checkPlot:
%         creates a new figure to check if the filtering worked
%     threshold; double [5]
%
%     removePixelAlerts: bool [true]
%         removes all pixels that were tagges as 'failed' (after version
%         2021.1.beta3)


arguments
   kwargs.nFiles = 'none';
   kwargs.checkPlot (1,1) {mustBeBoolean(kwargs.checkPlot)}= false
   kwargs.save {mustBeBoolean(kwargs.save)} = true;
   kwargs.removePixelAlerts = false;

   filterArgs.cutOff = 'none';
   filterArgs.includeHotPixel  = 0
   filterArgs.chi = 0
   filterArgs.winSize = 3
   filterArgs.threshold = 'none'; %in Gauss
end

% define default values for function
defaults = struct('threshold', 5);
filterArgs = ask_arguments(filterArgs, defaults);

% checks and detects if a path was given otherwise you can select one
nFiles = automatic_input_ui__(kwargs.nFiles, 'type', 'file', 'multiselect', 'off');

filteredMaps = {};

for i = 1:size(nFiles,2)
    iFile = nFiles{i};
    expData = load(iFile);
    msg = sprintf('Thresholding data for file << %s >>', iFile);
    logMsg('info',msg,1,0);
    
    if is_B111(expData)
        filterData = expData.B111ferro;
    else
        filterData = expData.Bz;
    end
    
    if kwargs.removePixelAlerts && sum(strcmp(fieldnames(expData), 'pixelAlerts'))
        sprintf('<>     removing << %i >> pixels that failed Lorentzian fits (see documentation)\n', numel(nonzeros(expData.pixelAlertss)));
        filterData = filterData .* ~expData.pixelAlerts;
    end
    
    if filterArgs.chi
        % check if cutOff is defined
        if strcmp(filterArgs.cutOff, 'none')
            msg = sprintf('''cutOff'' is set to ''none''. Threshold will not be defined according to the chi values.\n                       please set: e.g. ''cutOff'', 5');
            logMsg('warn',msg,1,0);
        end
        
        try
            chi = expData.chi2Neg1 + expData.chi2Neg2 + expData.chi2Pos1 + expData.chi2Pos2;
        catch
            msg = sprintf('expData does not contain chi2 values. Not using Chi2 filtering!');
            logMsg('error',msg,1,0);
        end
    else
        chi = 0;
    end
    
    % apply filter
    passFilterArgs = namedargs2cell(filterArgs);
    filterData = filter_hot_pixels(filterData, passFilterArgs{:});
    
	if is_B111(expData)
        expData.B111ferro_unfiltered = expData.B111ferro;
        expData.B111ferro = filterData;
    else
        expData.Bz_unfiltered = expData.Bz;
        expData.Bz = filterData;
    end
    
    filteredMaps{end+1} = expData;
    
    if kwargs.save
        % save data with new fileName
        suffix = sprintf('_thresh(%s).mat', num2str(filterArgs.threshold));
        iFileNew = strrep(iFile,'.mat', string(suffix));
        msg = sprintf('filtered data for file << %s >>', iFileNew);
        logMsg('saving',msg,1,0);
        save(iFileNew, '-struct', 'expData')
    end

    if kwargs.checkPlot
        threshold_checkPlot(filteredMaps{i})
    end
end
end

function threshold_checkPlot(data)
%threshold_checkPlot(data)
    fig = figure('Name', 'Threshold checkPlot', 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.3]);

    movegui(fig,'center')

    [~, dataName, ~] = is_B111(data);

    ax1 = subplot(1,3,1);
    QDM_figure(data.(sprintf('%s_unfiltered', dataName)), 'preThreshold', 'none', ...
        'ax', ax1, 'unit', 'muT','title','before Threshold');
    
    ax2 = subplot(1,3,2);
    QDM_figure(data.(dataName), 'preThreshold', 'none', ...
        'ax', ax2, 'unit', 'muT', 'title','after Threshold');
    
    ax3 = subplot(1,3,3);
    QDM_figure(data.(dataName) - data.(sprintf('%s_unfiltered', dataName)), ...
        'preThreshold', 'none', 'colormap', 'rwb', ...
        'ax', ax3, 'unit', 'muT', 'title','delta');
    linkaxes([ax1 ax2 ax3])
    
end