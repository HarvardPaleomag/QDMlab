function expData = make_square_hole(kwargs)
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
bData = expData.(dataName);

[bkgy, bkgx] = pick_box2('expData', expData, 'title', 'Select background value', 'point',true);
limx=max(1,bkgx-2);
limX=min(size(bData,2),bkgx+2);
limy=max(1,bkgy-2);
limY=min(size(bData,1),bkgy+2);
bkg=mean2(bData(limy:limY,limx:limX));

[row, col] = pick_box2('expData', expData, 'title', 'Select area to remove');

%cropping a B111 map set
% nReplace=0;
% for i = min(xx):max(xx)
%     for j = min(col):max(row)
%         bDataReplaced(j,i) = bkg;
%         nReplace = nReplace + 1;
%     end
% end
bData(row(1):row(2),col(1):col(2)) = bkg;

expData.([dataName '_original']) = expData.(dataName);
expData.(dataName) = bData;

if kwargs.cropFigure || kwargs.save
    fig = figure('Units', 'normalized', ...
                 'Position',[0.2 0.2 0.5 0.5], 'Name', 'replaced map');

    QDM_figure(bData, 'title', 'replaced map', 'fig', fig)
end

%% save data with new fileName
if kwargs.save
    [filepath,~,~]=fileparts(dataFile{1,1});
    
    iFileNew = strrep(dataFile{1,1}, '.mat','_Hole.mat');
    fprintf('<>     SAVING: replaced data to file << %s >>\n', iFileNew);
    
    saveas(fig,[filepath '/B111Cropped.png'])
    save(iFileNew,'-struct','expData');
end

end