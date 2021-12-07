function expData = replace_with_constant(kwargs)
%[expData] = replace_with_constant('filePath', 'save', 'editFigure')
%[expData] = subtract_constant('filePath', 'save', 'editFigure')
% This script takes an input Bz map, asks for a box, crops to that box, and
% outputs Bz and Bt maps, along with the accessory parameters

arguments
    kwargs.filePath = 'none'
    kwargs.save = true;
    kwargs.editFigure = true;
end

filePath = automatic_input_ui__(kwargs.filePath, 'type', 'file', 'title', 'Pick a magnetic field map file');
expData = load(filePath{:});

[~, dataName, ~] = is_B111(expData);
bData = expData.(dataName);

[nROI, coordinates] = pick_box(bData, 'title', 'select Area for replacement');
[lin, col] = pick_box2('expData', expData, 'title', 'Select constant to replace', 'point', true);

const = bData(lin,col);

for i = 1:size(nROI,2)
    coord = coordinates{i};
    x1 = coord(1);
    x2 = coord(1) + coord(3);
    y1 = coord(2);
    y2 = coord(2) + coord(4);

    bData(y1:y2, x1:x2) = const;
end

if kwargs.editFigure
    fig = figure('Units', 'normalized', ...
                 'Position',[0.2 0.2 0.5 0.5], 'Name', 'edited map');
    QDM_figure(bData, 'title', 'edited map', 'fig', fig);
    
    if kwargs.save
        [path,name] = fileparts(filePath);
        imagePath = fullfile(path, sprintf('%s_%sRep.png', name, dataName));
        msg = sprintf('saving image << %s >>', imagePath');
        logMsg('info',msg,1,0);
        saveas(fig, imagePath)
    end
end

if kwargs.save
    % save edited map to struct 
    expData.([dataName '_original']) = expData.(dataName);
    expData.(dataName) = bData;
    
    iFileNew = strrep(filePath{1,1}, '.mat', sprintf('_%sRep.mat', dataName));
    msg = sprintf('saving edited data to file << %s >>', iFileNew');
    logMsg('info',msg,1,0);
    save(iFileNew,'-struct','expData');
end

end
