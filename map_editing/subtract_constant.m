function expData = subtract_constant(kwargs)
%[expData] = subtract_constant('editFigure', 'filePath', 'save')
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

    [fPath,fileName,~]=fileparts(filePath{1,1});

    iFileNew = strrep(filePath{1,1}, '.mat', sprintf('_%sSub.mat', dataName));
    fprintf('<>     SAVING: cropped data to file << %s >>\n', iFileNew);
    
    saveas(fig,fullfile(fPath, sprintf('%s_%sSub.png', fileName, dataName)))
    save(iFileNew,'-struct','expData');
end

end
