function save_transformed_map(movingPath, referenceData, kwargs)
%save_transformed_map(movingPath, referenceData; 'reference', 'save', 'checkPlot')
% 
% Parameters
% ----------

arguments
    movingPath
    referenceData
    kwargs.reference = 'LED'
    kwargs.save = true
    kwargs.checkPlot = 1
end

movingData = load(movingPath);
newFile = movingData;

[~, dataName, ledName] = is_B111(movingData);

[transForm, refFrame] = get_image_tform(movingData.(kwargs.reference), ...
                        referenceData,...
                        'checkPlot', kwargs.checkPlot);
                    
 newFile.(ledName) = tform_data(movingData.(ledName), transForm, refFrame);
 newFile.(dataName) = tform_data(movingData.(dataName), transForm, refFrame);
 newFile.reference = kwargs.reference;
 newFile.referenceData = referenceData;
 
 if kwargs.save    
    iFileNew = strrep(movingPath, '.mat','_Aligned.mat');
    fprintf('<>     SAVING: cropped data to file << %s >>\n', iFileNew);
    save(iFileNew,'-struct','newFile');
end
 
 
