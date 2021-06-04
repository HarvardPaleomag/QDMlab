function [results, files, nROI, nMasks] = estimate_coercivity(nFolders, kwargs, selection, error, filter)
%[results, files, nROI, nMasks] = estimate_coercivity(nFolders; 'fileName', 'transFormFile', 'fixedIdx', 'upCont', 'checkPlot', 'reverse', 'nROI', 'freeHand', 'freeHandSelection', 'selectionThreshold', 'bootStrapN', 'pixelError', 'filterProps')
%DEPRECATED: please use demag_behavior

arguments
    nFolders cell {foldersMustExist(nFolders)}
    kwargs.fileName char {fileMustExistInFolder(kwargs.fileName, nFolders)} = 'Bz_uc0'
    kwargs.transFormFile = 'none'
    kwargs.fixedIdx (1,1) {mustBePositive} = 1
    kwargs.upCont = false;
    
    kwargs.checkPlot  (1,1) {mustBeBoolean(kwargs.checkPlot)} = false;
    kwargs.reverse  (1,1) {mustBeBoolean(kwargs.reverse)} = false;
    kwargs.nROI = false
    
    selection.freeHand  (1,1) {mustBeBoolean(selection.freeHand)} = false
    selection.freeHandSelection (1,1) {mustBeBoolean(selection.freeHandSelection)} = false
    selection.selectionThreshold (1,1) {mustBeNumeric} = 0.25
    
    error.bootStrapN = 1
    error.pixelError = 4
    
    filter.filterProps struct = struct();

end

msg = sprintf('please use demag_behavior. estimate_coercivity will be removed in a later update');
logMsg('deprecated',msg,1,0);

kwargs = namedargs2cell(kwargs);
selection = namedargs2cell(selection);
error = namedargs2cell(error);
filter = namedargs2cell(filter);

[results, files, nROI, nMasks] = demag_behavior(nFolders, kwargs{:}, selection{:}, error{:}, filter{:});
