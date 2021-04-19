function [row, col] = pick_box2(kwargs)
%
% Parameters
% ----------
%     expData: dataStruct or filepath ['none']
%       QDM data structure or the file path to the data. 
%       If 'none', ui lets you pickj the data file and loads it.
%     title: str, ['Pick Area (+/- to change colorscale)']
%       Title of the plot
%
% Returns
% -------
%   row, col indices of the box
arguments
    kwargs.expData = 'none';
    kwargs.title = 'Pick Area (+/- to change colorscale)';
end

%% get data
if strcmp(kwargs.expData, 'none')
    dataFile = automatic_input_ui__(kwargs.dataFile, 'type', 'file', 'title', 'Pick a magnetic field map file', 'single', true);
    expData = load(dataFile);
else
    expData = kwargs.expData;
end

[~, dataName, ~] = is_B111(expData);
bData = expData.(dataName);
% ledData = expData.(ledName);

%% Crop a single region for source fitting
fprintf('<>   select area:  (+) saturate color scale, (-) desaturate color scale,\n')
fprintf('<>                 (*) or (/) restore original color scale\n');

pickFigure = figure();
ax = QDM_figure(bData, 'fig', pickFigure, 'title', kwargs.title, 'return', 'ax');

hold(ax, 'on')

flag = 1;

%% pick region and +/- climits
row = [];
col = [];
while flag
    [c, l, key] = ginput(1);
    switch key
        case 43
            caxis(caxis/sqrt(10));
        case 45
            caxis(caxis*sqrt(10));
        case {42, 47}
            caxis auto, caxis([-1, 1]*max(abs(caxis)));
        case 1
            row = [row; round(l)];
            col = [col; round(c)];
            plot(ax, col, row, '+m','MarkerSize',12)
        case 27
            return
    end
    flag = length(row) ~= 2;
end

row = sort(row, 1);
col = sort(col, 1);

%% check indices
if row(1) <= 0
    row(1) = 1;
end
if col(1) <= 0
    col(1) = 1;
end

if row(2) > size(bData, 1)
    row(2) = size(bData, 1);
end
if col(2) > size(bData, 2)
    col(2) = size(bData, 2);
end

pause(1)
close(pickFigure)