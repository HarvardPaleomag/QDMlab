function expData = subtract_constant(kwargs)
%[expData] = subtract_constant('filePath', 'save', 'editFigure')
% This script takes an input Bz map, asks for a box, crops to that box, and
% outputs Bz and Bt maps, along with the accessory parameters
% 
% Parameters
% ----------
%   filePath: ('none')
%   save: (true)
%   editFigure: (true)
% 
% Returns
% ----------
%   expData:

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

msg = sprintf('Constant %s value of %.2e subtracted from map', dataName, const);
logMsg('info',msg,1,0);

[fPath,fileName,~]=fileparts(filePath{1,1});

if kwargs.editFigure
    fig = figure('Units', 'normalized', ...
                 'Position',[0.2 0.2 0.5 0.5], 'Name', 'edited map');
    QDM_figure(bData, 'title', 'edited map', 'fig', fig);
    
    if kwargs.save
        imagePath = fullfile(fPath, sprintf('%s_%sSub.png', fileName, dataName));
        msg = sprintf('saving image << %s >>', imagePath');
        logMsg('info',msg,1,0);
        saveas(fig, imagePath)
    end
end

if kwargs.save
    % save edited map to struct 
    expData.([dataName '_original']) = expData.(dataName);
    expData.(dataName) = bData;
    
    iFileNew = strrep(filePath{1,1}, '.mat', sprintf('_%sSub.mat', dataName));
    msg = sprintf('saving edited data to file << %s >>', iFileNew');
    logMsg('info',msg,1,0);
    save(iFileNew,'-struct','expData');
end

end
