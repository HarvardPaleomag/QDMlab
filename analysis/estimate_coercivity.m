function [results, files, nROI, nMasks] = estimate_coercivity(nFolders, kwargs, selection, error, filter)%demag_behavior
%[results, files, nROI, nMasks] = estimate_coercivity(nFolders; 'bootStrapN', 'checkPlot', 'fileName', 'filterStruct', 'fixedIdx', 'freeHand', 'freeHandSelection', 'nROI', 'pixelError', 'reverse', 'selectionThreshold', 'threshold', 'transFormFile', 'upCont')
% These codes (1) register the maps and (2) analizes a user selected magnetic
% pattern for changes from one map to the next.(folders, varargin)
%
% Parameters
% ----------
%     folders: list
%         list of absolute path to each data folder. First entry in the list
%         is used as reference image.
%         Note: calculations in the order the folders are given
%     fileName: str ['Bz_uc0.mat']
%         name of the .mat file to be used.
%     transFormFile: str ['none']
%         absolute path of a file where previously calculated transformations
%         are in. If file does not exist the file will be created so that
%         you dont have to do this every time.
%     removeHotPixel: bool [0]
%         if True hotpixels will be removed
%         if False hotpixels will not be removed
%     includeHotPixel: bool [0]
%         | if 1: hotpixel value is also used to calculate the new value that replaces the hot pixel
%         | if 0: only pixels in window with winSize are used to calculate the new value that replaces the hot pixel
%     selectionThreshold: numeric [0.25]
%         | defines the Threshold above which the mask is created.
%         | **Example:** selectionThreshold = 0.5
%         | :code:`maskSelection = [1 2 0; 1 1 2; 0 1 1]`  ->
%         | :code:`mask          = [0 1 0; 1 0 1; 0 0 1]`
%     upCont: cell
%         Cellarray of upward continuation distances for each of the files.
%         Needs to be in same order.
%     filter: str
%         Uses filter_hot_pixels to clean the maps before calculating the
%         results.
%         Note: 'includeHotPixel' is false by default
%     freeHand: bool [0]
%         Instead of using a rectangular selection use freehand drawn
%     freeHandSelection: bool [0]
%         | If true: the selectionThreshold is used to create the mask from the maskSelection.
%         | If false: the mask = maskSelection
%     fixedIdx: int [1]
%         index of the reference LED image. This will fixed while the other
%         image is transformed.
%     chi:
%         | if true the chi2 value is used to filter and not the data itself
%         | **Note:** the chi2 value is calculated from the sum(pos1, pos2, neg1,neg2) chi values.
%     nROI: cell
%         | cell with region of interest.
%         | **Note:** the ROI and the selection are not necessarily the same.|
%           If the ROI was selected with freeHand then ROI == selection.
%           Otherwise it is a rectangle around the selection.
%     reverse: bool ==> CURRENTLY NOT SUPPORTED
%         default: false
%         if true:  refernce - tform -> target
%         if false: target   - tform -> reference
%     bootStrapN: int, bool [1]
%         This uses a boot strapping approach to estimate errors. It will
%         shift the mask by 'pixelError' 'bootStrapN' times. The
%         result value is calculated from the mean of all calculations and an
%         error is estimated.
%         If 'bootStrapError' == 1 only one value is calculated -> no error
%         estimation.
%     pixelError: int [0]
%         Used for bootstrapping an error estimate.
%         The number of pixels the mask will be randomly shifted to estimate
%         an error.
%
% Returns
% -------
%     struct with 7 cells for each result (iFiles,
%     pPixels, pPixelRats, iPixelThreshs,
%     sumDiffs, errs, maskSums, iMasks, transData).
%     Each is a cell with result{i,j} with ith mask and jth file:
%
% hint
% ----
%     **Example:** Here we have a cell of folders, containing QDM images ('B111dataToPlot', 'Bz_uc0').
%     literal blocks::
%         >> nFolders = {'home/NRM, 'home/10mT','home/20mT','home/30mT'};
%
%     The most simple way to estimate the coercivity is to call:
%     literal blocks::
%         >> estimate_coercivity(nFolders)
%
%     This calls the function with the default parameters (the values in square brackets *[]*).     


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
