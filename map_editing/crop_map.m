function [expData, row, col] = crop_map(kwargs)
% This script takes an input Bz map, asks for a box, crops to that box, and
% outputs Bz and Bt maps, along with the accessory parameters

arguments
    kwargs.filePath = 'none';
    kwargs.save = true;
    kwargs.checkPlot = true;
    kwargs.row = 'none';
    kwargs.col = 'none';
end

if isstruct(kwargs.filePath)
    expData = kwargs.filePath;
    filePath = expData.filePath;
else
    filePath = automatic_input_ui__(kwargs.filePath, 'type', 'file', 'title', 'Pick a magnetic field map file', 'single', true);
    expData = load(filePath);
end

[~, dataName, ~] = is_B111(expData);

[row, col] = pick_box2('expData', expData, 'title', 'Select area to crop', 'even', true);

bData = expData.(dataName);
%cropping a B111 map set
corners=[row(1),row(2);col(1),col(2)];
bDataCropped=bData(row(1):row(2),col(1):col(2));

expData.(dataName) = bDataCropped;
expData.([dataName '_original']) = bData;
expData.corners = corners;

if kwargs.checkPlot || kwargs.save
    fig = figure('Units', 'normalized', ...
                 'Position',[0.2 0.2 0.5 0.5], 'Name', 'cropped map');

    QDM_figure(bDataCropped, 'title', 'cropped map', 'fig', fig)
end

%% save data with new fileName
if kwargs.save
    [filePath,fileName,~]=fileparts(filePath);
    
    iFileNew = strrep(filePath, '.mat','_Cropped.mat');
    fprintf('<>     SAVING: cropped data to file << %s >>\n', iFileNew);
    saveas(fig,fullfile(filePath, sprintf('%s_Hole.png', fileName)))
    save(iFileNew,'-struct','expData');
end

end