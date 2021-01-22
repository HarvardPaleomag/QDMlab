function filteredData = filter_hot_pixels(data, kwargs)
% This function takes a B111ferro file and filters it by replacing the hot 
% pixel the with the mean of the surrounding pixels (7x7). 
% %todo add size to arguments
% 
% Parameters
% ----------
%     positional
%     ----------
%     data:
%         data matrix to be filtered
%     
%     optional parameters
%     -------------------
%     cutOff: 
%         how many standard deviations have to be exceeded for the pixel to
%         be filtered.
%         optional, default = 4
%     includeHotPixel:
%         if true: the mean will be calculated including the hot pixel
%         if false: the mean is calculated after setting the pixel to nan
%     chi: array
%         if chi is provided the data will be filtered according to the chi
%         values.
%     winSize: int
%         specifies the number of pixels to the left AND right to be used for
%         averaging
%         if nan: values are replaced by nan
%     checkPlot:
%         creates a new figure to check if the filtering worked


arguments
   data
   kwargs.cutOff = 'none';
   kwargs.includeHotPixel = 0
   kwargs.checkPlot = 0
   kwargs.chi = 0
   kwargs.winSize = 3
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
aboveStd = abs(data) >= 5;
data(aboveStd) = nan;

%% chose mode for hot pixel calculation
if ~strcmp(cutOff, 'none')
    fprintf('<>   FILTER: data where data/chi is %i standard deviations above the median\n', cutOff)

    if chi ~= 0
        % pefilter chi values to catch extreme outlier
        chiFilter = abs(chi) > nanmean(chi, 'all') + 15 * nanstd(chi, 0, 'all');
        chi(chiFilter) = nan;

        % calculate the mean over all pixels
        dMed = nanmedian(abs(chi), 'all');

        %calculate the standard deviation of all pixels
        dStd = nanstd(chi, 0, 'all');
        disp(['<>           using chi2 values: median = ' num2str(dMed) '; std = ' num2str(dStd)])

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
    end
end

%% calculate pixel substitution
% copy original data into new array
filteredData = data;

% initiate filtered pixel array for plotting
filteredPixels = zeros(dshape);

dataMedian = nanmedian(abs(data), 'all');

% replace poixels with nan if specified
if isnan(winSize)
    % set pixel value to nan
    filteredData(aboveStd) = nan;
    filteredPixels(aboveStd) = 1;
% otherwise calculate the mean over winSize
else
    for row = 1:dshape(1)
        for col = 1:dshape(2)
            % pixels that exceed the mean + 4 std deviations are replaced by
            % the mean of a square of pixels 7x7 pixels
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
    fprintf('<>   FILTER: B > +- 5G: removed %i / %i pixel = %.2f%%\n', n_pixels, numel(filteredPixels), n_pixels/numel(filteredPixels)*100)
else
    fprintf('<>   FILTER: by %i stdev: removed %i / %i pixel = %.2f%% | median = %.2e, std = %.2e\n', cutOff, n_pixels, numel(filteredPixels), n_pixels/numel(filteredPixels)*100, dMed, dStd)
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
end