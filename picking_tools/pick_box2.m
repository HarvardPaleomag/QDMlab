function [lin, col] = pick_box2(kwargs)
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
lin = [];
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
            lin = [lin; round(l)];
            col = [col; round(c)];
            plot(ax, col, lin, '+m','MarkerSize',12)
        case 27
            return
    end
    flag = length(lin) ~= 2;
end

lin = sort(lin, 1);
col = sort(col, 1);

%% check indices
if lin(1) <= 0
    lin(1) = 1;
end
if col(1) <= 0
    col(1) = 1;
end

if lin(2) > size(bData, 1)
    lin(2) = size(bData, 1);
end
if col(2) > size(bData, 2)
    col(2) = size(bData, 2);
end

pause(1)
close(pickFigure)