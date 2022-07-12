function save_resized_raw_data(inFolder, outFolder, binSize, check)
%save_resized_raw_data(inFolder, outFolder, binSize; 'check')
% Saves newly binned data to disc.
%
% Parameters
% ----------
%   inFolder:
%   outFolder:
%   binSize:
%   check:
%   check: (True)

    arguments
        inFolder
        outFolder
        binSize
        check = True
    end

    if inFolder == outFolder
        msg = sprintf('inFolder and outFolder are the same !! Not allowed, because data would be overwritten.');
        logMsg('error',msg,1,0);
    end
    dataFiles = dir(fullfile(inFolder,'run_0000*.mat'));
    for i = 1: size(dataFiles, 1)
        dataFile = dataFiles(i);
        expData = load(fullfile(dataFile.folder, dataFile.name));
        dNew = resize_raw_data(expData, binSize, check);
        save(fullfile(outFolder, dataFile.name), '-struct', 'dNew')
    end

end

function dNew = resize_raw_data(expData, binSize, check)
%[dNew] = resize_raw_data(expData, binSize, check)
    arguments
        expData
        binSize
        check
    end
    %%
    expData1 = QDMreshape(expData.imgStack1, expData.imgNumRows, expData.imgNumCols);
    expData2 = QDMreshape(expData.imgStack2, expData.imgNumRows, expData.imgNumCols);
    %%
    eD1_resized = imresize(expData1, 1/binSize, 'method', 'box');
    eD2_resized = imresize(expData2, 1/binSize, 'method', 'box');
    %%
    eD1_resized_ = QDMreshape_reverse(eD1_resized, 51);
    eD2_resized_ = QDMreshape_reverse(eD2_resized, 51);
    %%
    if check
        eD1_resized_check = QDMreshape_reverse(expData1, 51);

        if all(expData.imgStack1 == eD1_resized_check, 'all')
            msg = sprintf('The new imagestack is transformed correctly');
            logMsg('debug',msg,1,0);
        else
            msg = sprintf('The new imagestack is transformed NOT correctly');
            logMsg('error',msg,1,0);
        end

    end
    %%
    dNew = expData;
    dNew.imgStack1 = eD1_resized_;
    dNew.imgStack2 = eD2_resized_;
    dNew.imgNumCols = size(eD1_resized, 2);
    dNew.imgNumRows = size(eD1_resized, 1);
    dNew = rmfield(dNew,'filePath');
    
    if check
        
        subplot(2,2,1)
        title('original data')
        imagesc(expData1(:,:,1))

        subplot(2,2,2)
        title('resized data')
        imagesc(eD1_resized(:,:,1))

        subplot(2,2,3)
        title('resized and reshaped manually')
        dReshapeReverse = QDMreshape(eD1_resized_, dNew.imgNumRows, dNew.imgNumCols);
        imagesc(dReshapeReverse(:,:,1))

        subplot(2,2,4)
        title('resized and reshaped prepare_raw_data')
        dReshapePrepare = prepare_raw_data(dNew,1,1);
        imagesc(dReshapePrepare(:,:,1))
    end

end