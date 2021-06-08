function filteredData = filter_hot_pixels(data, kwargs)
%[filteredData] = filter_hot_pixels(data; 'cutOff', 'includeHotPixel', 'checkPlot', 'chi', 'winSize', 'threshold')
<<<<<<< Updated upstream
% This function takes a B111ferro file and filters it by replacing the hot 
% pixel the with the mean of the surrounding pixels (7x7). 
% %todo add size to arguments
% 
% Parameters
% ----------
%     data: double
%         data matrix to be filtered
%     cutOff: double ['none']
%         how many standard deviations have to be exceeded for the pixel to
%         be filtered.
%     includeHotPixel: bool [false]
%         if true: the mean will be calculated including the hot pixel
%         if false: the mean is calculated after setting the pixel to nan
%     chi: array [0]
%         if chi is provided the data will be filtered according to the chi
%         values.
%     winSize: int [3]
%         specifies the number of pixels to the left AND right to be used for
%         averaging
%         if nan: values are replaced by nan
%         if 0: values are replaced by 0
%     checkPlot: bool [false]
%         creates a new figure to check if the filtering worked


arguments
   data
   kwargs.cutOff = 'none';
   kwargs.includeHotPixel {mustBeBoolean(kwargs.includeHotPixel)} = 0
   kwargs.checkPlot (1,1) {mustBeBoolean(kwargs.checkPlot)}= false
   kwargs.chi = 0
   kwargs.winSize = 3
   kwargs.threshold = 5
end

% define optional arguments
cutOff = kwargs.cutOff;
includeHotPixel = kwargs.includeHotPixel;
checkPlot = kwargs.checkPlot;
chi = kwargs.chi;
winSize = kwargs.winSize;

% shape of the data (row,col)
dshape = size(data);

% pefilter data values to catch extreme outlier
if ~isequal(kwargs.threshold, false)
    aboveStd = abs(data) >= kwargs.threshold;
    data(aboveStd) = nan;
else
    aboveStd = zeros(size(data));
end

if kwargs.threshold && numel(nonzeros(aboveStd))
    msg = sprintf('found %i pixels with field values above %.1e G', numel(nonzeros(aboveStd)), kwargs.threshold);
    logMsg('debug',msg,1,0);
end




%% Check if cutOff value is  given
if strcmp(cutOff, 'none') && ~isequal(chi, false)
    msg = sprintf('filtering according to Chi alues, but no ''cutOff'' value set. Please check add: e.g. << ''cutOff'', 5 >> !');
    logMsg('error',msg,1,0);
end

%% chose mode for hot pixel calculation
if ~strcmp(cutOff, 'none')
    if all(size(chi) == dshape) 
        dtype = 'chi';
    else
        dtype = 'data';
    end
    msg = sprintf('filtering where %s is %i standard deviations above the median', dtype, cutOff);
    logMsg('info',msg,1,0);

    if all(size(chi) == dshape)
        % pefilter chi values to catch extreme outlier
        chiFilter = abs(chi) > nanmean(chi, 'all') + 15 * nanstd(chi, 0, 'all');
        chi(chiFilter) = nan;

        % calculate the mean over all pixels
        dMed = nanmedian(abs(chi), 'all');

        %calculate the standard deviation of all pixels
        dStd = nanstd(chi, 0, 'all');
        
        msg = sprintf('using chi2 values: median = %.2e; std = %.2e',dMed, dStd);
        logMsg('debug',msg,1,0);

        %create boolean array with pixels with intensity higher than cutoff
        aboveStd = aboveStd | (chi > dMed + cutOff * dStd);
        aboveStd = aboveStd | chiFilter;
    else
        % calculate the mean over all pixels
        dMed = nanmedian(data, 'all');

        %calculate the standard deviation of all pixels
        dStd = nanstd(data, 0, 'all');
        %create boolean array with pixels with intensity higher than cutoff
        aboveStd = aboveStd | abs(data) > dMed + cutOff * dStd;
        
        msg = sprintf('using data values: median = %.2e; std = %.2e',dMed, dStd);
        logMsg('debug',msg,1,0);
    end
end

%% calculate pixel substitution
% copy original data into new array
filteredData = data;

% initiate filtered pixel array for plotting
filteredPixels = zeros(dshape);

dataMedian = median(abs(data), 'all', 'omitnan');

% replace pixels with nan if specified
if isnan(winSize)
    % set pixel value to nan
    filteredData(aboveStd) = nan;
    filteredPixels(aboveStd) = 1;
% replace pixels with 0 if specified
elseif winSize == 0
    % set pixel value to nan
    filteredData(aboveStd) = 0;
    filteredPixels(aboveStd) = 1;
% otherwise calculate the mean over winSize
else
    for row = 1:dshape(1)
        for col = 1:dshape(2)
            % pixels that exceed the mean + cutOff std deviations are replaced by
            % the mean of a square of winSize pixels
            if aboveStd(row, col)
                % set pixel to 1 to check which pixel were removed
                filteredPixels(row, col) = 1;

                if includeHotPixel == false
                    % set pixel value to nan
                    filteredData(row,col) = nan;
                end

                % set maximum and minimum of average square
                minrow = row-winSize;
                maxrow = row+winSize;

                mincol = col-winSize;
                maxcol = col+winSize;

                % deal with edge cases by reducing the window
                if minrow < 1
                    minrow = 1;
                end

                if mincol < 1 
                    mincol = 1;
                end


                if maxrow > dshape(1)
                    maxrow = dshape(1);
                end

                if maxcol > dshape(2)
                    maxcol = dshape(2);
                end

                % get the average window
                window = filteredData(minrow:maxrow, mincol:maxcol);
            
                % calculate the mean 
                % if include_hot_pixel without using the pixel itself
                new_val = nanmean(window, 'all');

                filteredData(row,col) = new_val;

            end
        end
    end
end

n_pixels = sum(sum(filteredPixels));

if strcmp(cutOff, 'none')
    if n_pixels
        msg = sprintf('B > +- %.3fG: removed %i / %i pixel = %.0f precent', kwargs.threshold, n_pixels, numel(filteredPixels), n_pixels/numel(filteredPixels)*100');
        logMsg('info',msg,1,0);
    end
else
    if n_pixels
        msg = sprintf('filtered by %i stdev: removed %i / %i pixel = %.2f%% | median = %.2e, std = %.2e', cutOff, n_pixels, numel(filteredPixels), n_pixels/numel(filteredPixels)*100, dMed, dStd');
        logMsg('info',msg,1,0);
    end
end

%% checkplot
if checkPlot
    figure
    % add hot pixel plot
    imagesc(filteredPixels > 0 );
    caxis([-1 1]*max(abs(caxis)));
    axis xy, axis equal, axis tight, axis off
    set(gca,'Fontsize',14);
    title('hot pixels','Fontsize',14);
    colormap(jet);
=======
% This function takes a B111ferro file and filters it by replacing the hot 
% pixel the with the mean of the surrounding pixels (7x7). 
% %todo add size to arguments
% 
% Parameters
% ----------
%     data: double
%         data matrix to be filtered
%     cutOff: double ['none']
%         how many standard deviations have to be exceeded for the pixel to
%         be filtered.
%     includeHotPixel: bool [false]
%         if true: the mean will be calculated including the hot pixel
%         if false: the mean is calculated after setting the pixel to nan
%     chi: array [0]
%         if chi is provided the data will be filtered according to the chi
%         values.
%     winSize: int [3]
%         specifies the number of pixels to the left AND right to be used for
%         averaging
%         if nan: values are replaced by nan
%         if 0: values are replaced by 0
%     checkPlot: bool [false]
%         creates a new figure to check if the filtering worked


arguments
   data
   kwargs.cutOff = 'none';
   kwargs.includeHotPixel {mustBeBoolean(kwargs.includeHotPixel)} = 0
   kwargs.checkPlot (1,1) {mustBeBoolean(kwargs.checkPlot)}= false
   kwargs.chi = 0
   kwargs.winSize = 3
   kwargs.threshold = 5
end

% define optional arguments
cutOff = kwargs.cutOff;
includeHotPixel = kwargs.includeHotPixel;
checkPlot = kwargs.checkPlot;
chi = kwargs.chi;
winSize = kwargs.winSize;

% shape of the data (row,col)
dshape = size(data);

% pefilter data values to catch extreme outlier
if ~isequal(kwargs.threshold, false)
    toFilter = abs(data) >= kwargs.threshold;
else
    toFilter = zeros(size(data));
end

if kwargs.threshold && numel(nonzeros(toFilter))
    msg = sprintf('found %i pixels with field values above %.1e G', numel(nonzeros(toFilter)), kwargs.threshold);
    logMsg('debug',msg,1,0);
end

if any(toFilter, 'all')
    data(toFilter) = nan;
end


%% Check if cutOff value is  given
if strcmp(cutOff, 'none') && ~isequal(chi, false)
    msg = sprintf('filtering according to Chi alues, but no ''cutOff'' value set. Please check add: e.g. << ''cutOff'', 5 >> !');
    logMsg('error',msg,1,0);
end

%% chose mode for hot pixel calculation
if ~strcmp(cutOff, 'none')
    if all(size(chi) == dshape) 
        dtype = 'chi';
    else
        dtype = 'data';
    end
    msg = sprintf('filtering where %s is %i standard deviations above the median', dtype, cutOff);
    logMsg('info',msg,1,0);

    if all(size(chi) == dshape)
        % pefilter chi values to catch extreme outlier
        chiFilter = abs(chi) > nanmean(chi, 'all') + 15 * nanstd(chi, 0, 'all');
        chi(chiFilter) = nan;

        % calculate the mean over all pixels
        dMed = nanmedian(abs(chi), 'all');

        %calculate the standard deviation of all pixels
        dStd = nanstd(chi, 0, 'all');
        
        msg = sprintf('using chi2 values: median = %.2e; std = %.2e',dMed, dStd);
        logMsg('debug',msg,1,0);

        %create boolean array with pixels with intensity higher than cutoff
        toFilter = toFilter | (chi > dMed + cutOff * dStd);
        toFilter = toFilter | chiFilter;
    else
        % calculate the mean over all pixels
        dMed = nanmedian(data, 'all');

        %calculate the standard deviation of all pixels
        dStd = nanstd(data, 0, 'all');
        %create boolean array with pixels with intensity higher than cutoff
        toFilter = toFilter | abs(data) > dMed + cutOff * dStd;
        
        msg = sprintf('using data values: median = %.2e; std = %.2e',dMed, dStd);
        logMsg('debug',msg,1,0);
    end
end

%% calculate pixel substitution
% copy original data into new array
filteredData = data;

% initiate filtered pixel array for plotting
filteredPixels = zeros(dshape);

dataMedian = median(abs(data), 'all', 'omitnan');

% replace pixels with nan if specified
if isnan(winSize)
    % set pixel value to nan
    filteredData(toFilter) = nan;
    filteredPixels(toFilter) = 1;
% replace pixels with 0 if specified
elseif winSize == 0
    % set pixel value to nan
    filteredData(toFilter) = 0;
    filteredPixels(toFilter) = 1;
% otherwise calculate the mean over winSize
else
    for row = 1:dshape(1)
        for col = 1:dshape(2)
            % pixels that exceed the mean + cutOff std deviations are replaced by
            % the mean of a square of winSize pixels
            if toFilter(row, col)
                % set pixel to 1 to check which pixel were removed
                filteredPixels(row, col) = 1;

                if includeHotPixel == false
                    % set pixel value to nan
                    filteredData(row,col) = nan;
                end

                % set maximum and minimum of average square
                minrow = row-winSize;
                maxrow = row+winSize;

                mincol = col-winSize;
                maxcol = col+winSize;

                % deal with edge cases by reducing the window
                if minrow < 1
                    minrow = 1;
                end

                if mincol < 1 
                    mincol = 1;
                end


                if maxrow > dshape(1)
                    maxrow = dshape(1);
                end

                if maxcol > dshape(2)
                    maxcol = dshape(2);
                end

                % get the average window
                window = filteredData(minrow:maxrow, mincol:maxcol);
            
                % calculate the mean 
                % if include_hot_pixel without using the pixel itself
                new_val = nanmean(window, 'all');

                filteredData(row,col) = new_val;

            end
        end
    end
end

n_pixels = sum(sum(filteredPixels));

if strcmp(cutOff, 'none')
    if n_pixels
        msg = sprintf('B > +- %.1fG: removed %i / %i pixel = %.2f precent', kwargs.threshold, n_pixels, numel(filteredPixels), n_pixels/numel(filteredPixels)*100');
        logMsg('info',msg,1,0);
    end
else
    if n_pixels
        msg = sprintf('filtered by %i stdev: removed %i / %i pixel = %.2f%% | median = %.2e, std = %.2e', cutOff, n_pixels, numel(filteredPixels), n_pixels/numel(filteredPixels)*100, dMed, dStd');
        logMsg('info',msg,1,0);
    end
end

%% checkplot
if checkPlot
    figure
    % add hot pixel plot
    imagesc(filteredPixels > 0 );
    caxis([-1 1]*max(abs(caxis)));
    axis xy, axis equal, axis tight, axis off
    set(gca,'Fontsize',14);
    title('hot pixels','Fontsize',14);
    colormap(jet);
>>>>>>> Stashed changes
end