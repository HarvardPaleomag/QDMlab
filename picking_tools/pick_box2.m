function [row, col] = pick_box2(kwargs)
%[row, col] = pick_box2(<even>, <expData>, <point>, <title>)
%
% Parameters
% ----------
%     expData: dataStruct or filepath ['none']
%       QDM data structure or the file path to the data. 
%       If 'none', ui lets you pickj the data file and loads it.
%     title: str, ['Pick Area (+/- to change colorscale)']
%       Title of the plot
%     point: bool [false]
%       if true: you only pick one point instead of two
%     even: bool [false]
%       if true: returns row/col with even mod(2) dimensions
%
% Returns
% -------
%   row, col indices of the box

arguments
    kwargs.expData = 'none';
    kwargs.title = 'Pick Area (+/- to change colorscale)';
    kwargs.point {mustBeBoolean(kwargs.point)} = false;
    kwargs.even {mustBeBoolean(kwargs.even)} = false;
end

%% get data
if strcmp(kwargs.expData, 'none')
    filePath = automatic_input_ui__(kwargs.expData, 'type', 'file', 'title', 'Pick a magnetic field map file', 'single', true);
    expData = load(filePath);
else
    expData = kwargs.expData;
end

[~, dataName, ~] = is_B111(expData);
bData = expData.(dataName);
% ledData = expData.(ledName);

%% Crop a single region for source fitting
msg = sprintf('area:  (+) saturate color scale, (-) desaturate color scale');
logMsg('SELECT',msg,1,0);
msg = sprintf('       (*) or (/) restore original color scale');
logMsg('SELECT',msg,1,0);

pickFigure = figure('Units', 'normalized', ...
               'Position',[0.1 0.2 0.8 0.8],'NumberTitle', 'off', 'Name', 'Pick Point(s)');
           
[pickFigure, ax] = QDM_figure(bData, 'fig', pickFigure, 'title', kwargs.title);

hold(ax, 'on')

flag = 1;

%% pick region and +/- climits
row = [];
col = [];

while flag
    [c, r, key] = ginput(1); %used to be c,l 
    switch key
        case 43
            caxis(caxis/sqrt(10));
        case 45
            caxis(caxis*sqrt(10));
        case {42, 47}
            caxis auto, caxis([-1, 1]*max(abs(caxis)));
        case 1
            if kwargs.point
                col = round(c);
                row = round(r);
                flag = 0;
            else
                row = [row; round(r)];
                col = [col; round(c)];
                flag = length(row) ~= 2;
            end
            plot(ax, col, row, '+m','MarkerSize',12)
        case 27
            return
    end
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

if ~kwargs.point
    if row(2) > size(bData, 1)
        row(2) = size(bData, 1);
    end
    if col(2) > size(bData, 2)
        col(2) = size(bData, 2);
    end
end

if kwargs.even
    row = makeEven(row);
    col = makeEven(col);
end

pause(1)
close(pickFigure)

if kwargs.point
    msg = sprintf('returning point indices: x: %i; y: %i', col, row);
else
    msg = sprintf('returning box indices: lower left: (%i, %i); upper right: (%i, %i)', col(1), row(1),col(2), row(2));
end

logMsg('debug',msg,1,0);

end

function vals = makeEven(vals)
%[vals] = makeEven(vals)
    %  return if not 2 elements in vals
    if size(vals,2) ~= 2
        return
    end
    
    %enforce that the cropped array has even dimensions
    range=diff(vals);
    
    if ~mod(range,2)
        vals(2,1)=vals(2,1)-1;
    end
end