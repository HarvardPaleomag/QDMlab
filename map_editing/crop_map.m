function expData = crop_map(kwargs)
% This script takes an input Bz map, asks for a box, crops to that box, and
% outputs Bz and Bt maps, along with the accessory parameters

arguments
    kwargs.dataFile = 'none'
    kwargs.save = true;
    kwargs.cropFigure = true;
end

dataFile = automatic_input_ui__(kwargs.dataFile, 'type', 'file', 'title', 'Pick a magnetic field map file');
expData = load(dataFile{:});
[~, dataName, ~] = is_B111(expData);

[row, col] = pick_box2('expData', expData, 'title', 'Select area to crop');

%enforce that the cropped array has even dimensions
xrange=col(2)-col(1);
yrange=row(2)-row(1);
if mod(xrange,2)
    col(2,1)=col(2,1);
else
    col(2,1)=col(2,1)-1;
end
if mod(yrange,2)
    row(2,1)=row(2,1);
else
    row(2,1)=row(2,1)-1;
end


bData = expData.(dataName);
%cropping a B111 map set
corners=[row(1),row(2);col(1),col(2)];
bDataCropped=bData(row(1):row(2),col(1):col(2));

expData.(dataName) = bDataCropped;
expData.([dataName '_original']) = bData;

if kwargs.cropFigure || kwargs.save
    fig = figure('Units', 'normalized', ...
                 'Position',[0.2 0.2 0.5 0.5], 'Name', 'cropped map');

    QDM_figure(bDataCropped, 'title', 'cropped map', 'fig', fig)
end

%% save data with new fileName
if kwargs.save
    [filePath,fileName,~]=fileparts(dataFile{1,1});
    
    iFileNew = strrep(dataFile{1,1}, '.mat','_Cropped.mat');
    fprintf('<>     SAVING: cropped data to file << %s >>\n', iFileNew);
    saveas(fig,fullfile(filePath, sprintf('%s_Hole.png', fileName)))
    save(iFileNew,'-struct','expData');
end

end