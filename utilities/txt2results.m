function results = txt2results(folders,sourceID,kwargs)

arguments
    folders;
    sourceID;
    
    % general kwargs
    kwargs.fileName = false;
    kwargs.UC = 0;
    kwargs.tFormFile = false;
end

% format string for DipoleInversions.txt
% fstr = '%s %s %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f';
fstr = '%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f%f\t%f\t%f\t%f\t%f\t%f';

% get number of folders 
numberoffolders = length(folders);

% preallocate the cells:
preAllocatedArray = zeros(1,numberoffolders);
iFiles = cell(size(preAllocatedArray));
UC = preAllocatedArray;
moments = preAllocatedArray;
inclinations = preAllocatedArray;
declinations = preAllocatedArray;
heights = preAllocatedArray; 
dipolarity = preAllocatedArray;
xMin = preAllocatedArray;
xMax = preAllocatedArray;
yMin = preAllocatedArray;
yMax = preAllocatedArray;
xloc = preAllocatedArray;
yloc = preAllocatedArray;
fileResults = num2cell(preAllocatedArray);
tForms = cell(numberoffolders);
refFrames = cell(numberoffolders);

% crawl through file paths 
for i = 1:length(folders)
    folder = folders(i);
    
    
    % open DipoleInversions file
    fID = fopen(fullfile(folder{1},'DipoleInversions.txt'),'r');
    
    % find row in DipoleInversions file
    fgetl(fID);

    filetxt=textscan(fID,fstr);
    
    
    % find row corresponding to this sourceID
    row = find(contains(filetxt{1},sourceID));
    
    % make sure sourceID exists and there are no duplicates
    if isempty(row)
        error('Input:sourceNotFound',...
            'Error: supplied sourceID not found in %s',folder{1});
    elseif length(row) > 1
        error('Input:duplicateSourceID',...
            'Error: multiple fits found with same sourceID in %s',folder{1});    
    end
    
    
    if kwargs.fileName
        fileName = kwargs.fileName;
    else
        fileName = filetxt{2}{1};
    end
    
    % get filename of Bz map
    if contains(fileName, '.mat')
        iFile = fullfile(folder, filesep, fileName);
    else
        iFile = fullfile(folder, filesep, [fileName, '.mat']);
    end

    UC(i) = filetxt{3};
    moments(i) = filetxt{4};
    declinations(i) = filetxt{5};
    inclinations(i) = filetxt{6};
    heights(i) = filetxt{7};
    dipolarity(i) = filetxt{8};
    xloc(i) = filetxt{9};
    yloc(i) = filetxt{10};
    xMin(i) = min([filetxt{12},filetxt{11}]);
    xMax(i) = max([filetxt{12},filetxt{11}]);
    yMin(i) = min([filetxt{13},filetxt{14}]);
    yMax(i) = max([filetxt{13},filetxt{14}]);
    
    iFiles{i} = fullfile(folder{1},fileName);
    
    % add tForm and refFrame if desired
    if kwargs.tFormFile
        tForms{i} = nTransForms(iFile);
        refFrames{i} = nRefFrames(iFile);
    end
end

% build final results structure
results = struct();
results.nFiles = iFiles;

results.moments = moments;
results.decs = declinations;
results.incs = inclinations;
results.heights = heights;
results.dipolarity = dipolarity;
results.res = 1-dipolarity;
results.x = xloc;
results.y = yloc;
results.xMin = xMin;
results.xMax = xMax;
results.yMin = yMin;
results.yMax = yMax;


if kwargs.tFormFile
    results.tforms = tForms;
    results.refFrame = refFrames;
end

if kwargs.UC
    results.UC = kwargs.UC;
else
    tmp = strsplit(fileName,{'_uc','.'});
    
    mapUC = 0;
    for j = tmp
        if ~isempty(str2num(j{1}))
            mapUC = mapUC + str2num(j{1});
        end
    end
    
    %check if inferred UC looks good
    UCgood=input(['UC = ', num2str(mapUC+UC(1)),'? [Y]:'],'s');
    if isempty(UCgood)
        results.UC = mapUC + UC(1);
    elseif UCgood=='y' | UCgood=='Y'
        results.UC = mapUC + UC(1);
    else
        UCinput = str2num(input('Correct UC value (microns): ','s'));
        if isempty(UCinput)
            error('TypeError:UC',...
            'Error: UC not convertible to double');
        else
            results.UC = UCinput;
        end
        
    end

end

msg = sprintf('returning structure where result{i,j} is: ith ROI, jth file');
logMsg('result',msg,1,0);
