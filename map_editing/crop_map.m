function [expData, row, col] = crop_map(kwargs)
%[expData, row, col] = crop_map('filePath', 'save', 'checkPlot', 'row', 'col', 'title')
% This script takes an input Bz map, asks for a box, crops to that box, and
% outputs Bz and Bt maps, along with the accessory parameters

arguments
    kwargs.filePath = 'none';
    kwargs.save = true;
    kwargs.checkPlot (1,1) {mustBeBoolean(kwargs.checkPlot)}= false
    kwargs.row = false;
    kwargs.col = false;
    kwargs.even = true; %enforce even dimensions in crop to be compatible with downstream fft functions
    kwargs.title = 'Select area to crop';
end

if isstruct(kwargs.filePath)
    expData = kwargs.filePath;
    filePath = expData.filePath;
else
    filePath = automatic_input_ui__(kwargs.filePath, 'type', 'file', 'title', 'Pick a magnetic field map file', 'single', true);
    expData = load(filePath);
end

if ~kwargs.row & ~kwargs.col
    [row, col] = pick_box2('expData', expData, 'title', kwargs.title, 'even', true);
else
    row = kwargs.row;
    col = kwargs.col;
end

if kwargs.even
    if ~mod(row(2)-row(1),2)
        row(2)=row(2)-1;
    end
    if ~mod(col(2)-col(1),2)
        col(2)=col(2)-1;
    end
end

[~, dataName, ledName] = is_B111(expData);
bData = expData.(dataName);
led = expData.(ledName);

binning = detect_binning(expData);

%cropping a B map set
corners=[row(1),row(2);col(1),col(2)];
bDataCropped=bData(row(1):row(2),col(1):col(2));
ledCropped = led(row(1)*binning:row(2)*binning,col(1)*binning:col(2)*binning);

expData.(dataName) = bDataCropped;
expData.([dataName '_original']) = bData;
expData.(ledName) = ledCropped;
expData.([ledName '_original']) = led;
expData.corners = corners;

if kwargs.checkPlot || kwargs.save
    fig = figure('Units', 'normalized', ...
                 'Position',[0.2 0.2 0.5 0.5], 'Name', 'cropped map');
    QDM_figure(bDataCropped, 'title', 'cropped map', 'fig', fig)
end

%% save data with new fileName
if kwargs.save
    fullpath=filePath;
    [filePath,fileName,~]=fileparts(filePath);
    
    iFileNew = strrep(fullpath,'.mat','_Cropped.mat');
    fprintf('<>     SAVING: cropped data to file << %s >>\n', iFileNew);
    saveas(fig,fullfile(filePath, sprintf('%s_crop.png', fileName)))
    save(iFileNew,'-struct','expData');
end

end