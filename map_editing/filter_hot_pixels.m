function filteredData = filter_hot_pixels(data, varargin)
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

inParse = inputParser;
addRequired(inParse, 'data');
addParameter(inParse, 'cutOff', 3, @isnumeric);
addParameter(inParse, 'includeHotPixel', false, @islogical);
addParameter(inParse, 'checkPlot', false, @islogical);
addParameter(inParse, 'chi', false, @ismatrix);
addParameter(inParse, 'winSize', 3, @isnumeric);
parse(inParse, data, varargin{:});

% define optional arguments
cutOff = inParse.Results.cutOff;
includeHotPixel = inParse.Results.includeHotPixel;
checkPlot = inParse.Results.checkPlot;
chi = inParse.Results.chi;
winSize = inParse.Results.winSize;

% shape of the data (row,col)
dshape = size(data);

% pefilter data values to catch extreme outlier
aboveStd = abs(data) > nanmean(data, 'all') + 10 * nanstd(data, 0, 'all');
data(aboveStd) = nan;

%% chose mode for hot pixel calculation
if chi ~= 0
    % pefilter chi values to catch extreme outlier
    chiFilter = abs(chi) > nanmean(chi, 'all') + 10 * nanstd(chi, 0, 'all');
    chi(chiFilter) = nan;
    
    % calculate the mean over all pixels
    chiMean = nanmedian(abs(chi), 'all');

    %calculate the standard deviation of all pixels
    chiStd = nanstd(chi, 0, 'all');
    disp(['<>   using chi2 values: median = ' num2str(chiMean) '; std = ' num2str(chiStd)])

    %create boolean array with pixels with intensity higher than cutoff
    aboveStd = aboveStd | (chi > chiMean + cutOff * chiStd);
    aboveStd = aboveStd | chiFilter;
else
    % calculate the mean over all pixels
    dmean = nanmedian(data, 'all');

    %calculate the standard deviation of all pixels
    dstd = nanstd(data, 0, 'all');
    disp(['<>   using data values: median = ' num2str(dmean) '; std = ' num2str(dstd)])

    %create boolean array with pixels with intensity higher than cutoff
    aboveStd = aboveStd | abs(data) > dmean + cutOff * dstd;
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

%                 if abs(new_val) > 10 * dataMedian
%                     new_val = nan;
%                 end


                filteredData(row,col) = new_val;

            end
        end
    end
end

n_pixels = sum(sum(filteredPixels));
disp(['<>   filtering by ' num2str(cutOff) ' stdev: removed ' int2str(n_pixels) ' of ' int2str(numel(filteredPixels)) '. Corresponding to ', num2str((n_pixels/numel(filteredPixels))*100, 2) '%.'])

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