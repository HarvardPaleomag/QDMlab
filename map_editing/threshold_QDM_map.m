function threshold_QDM_map(nFiles, kwargs)
% threshold_QDM_map(nFiles, kwargs) loads a file, removes all pixels above
% a threshold and saves it with '_thresh' tag.
% Parameters
% ----------
%     nFolders: cell
%     cutOff: int [4]
%         how many standard deviations have to be exceeded for the pixel to
%         be filtered.
%     includeHotPixel: bool [false]
%         if true: the mean will be calculated including the hot pixel
%         if false: the mean is calculated after setting the pixel to nan
%     chi: array [false]
%         if chi is provided the data will be filtered according to the chi
%         values.
%     winSize: int [3]
%         specifies the number of pixels to the left AND right to be used for
%         averaging
%         if nan: values are replaced by nan
%     checkPlot:
%         creates a new figure to check if the filtering worked
%     remove_failed_pixels: bool [true]
%         removes all pixels that were tagges as 'failed' (after version
%         2021.1.beta3)


arguments
   nFiles;
   kwargs.cutOff = 'none';
   kwargs.includeHotPixel  = 0
   kwargs.checkPlot = 0
   kwargs.chi = 0
   kwargs.winSize = 3
   kwargs.remove_failed_pixels = true;
end

nFiles = correct_cell_shape(nFiles);

for i = 1:size(nFiles,2)
    iFile = nFiles{i};
    expData = load(iFile);
    fprintf('<>   Thresholding data for file << %s >>\n', iFile);

    if is_B111(expData)
        filterData = expData.B111ferro;
    else
        filterData = expData.Bz;
    end
    
    if kwargs.remove_failed_pixels && sum(strcmp(fieldnames(expData), 'failed_pixel'))
        sprintf('<>     removing << %i >> pixels that failed Lorentzian fits (see documentation)\n', numel(nonzeros(expData.failedPixels)));
        filterData = filterData .* ~expData.failedPixels;
    end
    
    if kwargs.chi
        chi = expData.chi2Neg1 + expData.chi2Neg2 + expData.chi2Pos1 + expData.chi2Pos2;
    else
        chi =0;
    end
    
    filterData = filter_hot_pixels(filterData, 'cutOff',kwargs.cutOff, ...
                                     'includeHotPixel', kwargs.includeHotPixel,...
                                     'chi',chi,'winSize', kwargs.winSize);
    
    suffix = sprintf('_thresh(%s).mat', num2str(kwargs.cutOff));
    iFileNew = strrep(iFile,'.mat', string(suffix));
    fprintf('<>     SAVING: filtered data for file << %s >>\n', iFileNew);
end