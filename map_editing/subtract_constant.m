function expData = subtract_constant(kwargs)
% This script takes an input Bz map, asks for a box, crops to that box, and
% outputs Bz and Bt maps, along with the accessory parameters

arguments
    kwargs.dataFile = 'none'
    kwargs.save = true;
    kwargs.editFigure = true;
end

dataFile = automatic_input_ui__(kwargs.dataFile, 'type', 'file', 'title', 'Pick a magnetic field map file');
expData = load(dataFile{:});

[~, dataName, ~] = is_B111(expData);
bData = expData.(dataName);

[lin, col] = pick_box2('expData', expData, 'title', 'Select constant to subtract', 'point', true);

const = bData(lin,col);
bData = bData - const;

fprintf('<>   INFO: Constant %s value of %.2f subtracted from map', dataName, const);

if kwargs.editFigure
    fig = figure('Units', 'normalized', ...
                 'Position',[0.2 0.2 0.5 0.5], 'Name', 'edited map');
    QDM_figure(bData, 'title', 'cropped map', 'fig', fig)
end

if kwargs.save
    % save edited map to struct 
    expData.([dataName '_original']) = expData.(dataName);
    expData.(dataName) = bData;

    [filePath,fileName,~]=fileparts(dataFile{1,1});

    iFileNew = strrep(dataFile{1,1}, '.mat', sprintf('_%sSub.mat', dataName));
    fprintf('<>     SAVING: cropped data to file << %s >>\n', iFileNew);
    
    saveas(fig,fullfile(filePath, sprintf('%s_%sSub.png', fileName, dataName)))
    save(iFileNew,'-struct','expData');
end

end
