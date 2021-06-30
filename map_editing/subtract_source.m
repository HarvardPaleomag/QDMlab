function expData = subtract_source(kwargs)
%[expData] = subtract_source('filePath', 'save', 'fitOrder', 'checkPlot')
% This script takes an input Bz map, asks for a box, crops to that box, and
% outputs Bz and Bt maps, along with the accessory parameters
%
% Parameters
% ----------
%     filePath: str ['none']
%     save: bool [true]
%     fitOrder: 1, 2 or 3 [1]
%       1=dipole; 2=quadrupole; 3=octapole
%
% Returns
% -------
%     Map: double
%       after subtraction
%% check why weird
arguments
    kwargs.filePath = 'none'
    kwargs.save = true;
    kwargs.fitOrder = 1;
    kwargs.checkPlot (1,1) {mustBeBoolean(kwargs.checkPlot)}= false
end

show_references()

filePath = automatic_input_ui__(kwargs.filePath, 'type', 'file', ...
    'title', 'Pick a magnetic field map file', 'single', true);
expData = load(filePath);

[path, fileName, ext] = fileparts(filePath);
expData.filePath = path;

[b111, dataName, ~] = is_B111(expData);

if b111
    error('B111 file selected, please use Bz file')
end

bData = expData.(dataName);
expData.([dataName, '_original']) = bData;

[~, row, col] = crop_map('filePath', expData, 'save', false, 'checkPlot', false, 'title', 'select source');
%% maps ending in _dip for naming the file
dipMaps = dir(sprintf('/Users/mike/Dropbox/science/_projects/QDMlab_paper/data/NRM/4x4Binned/*-dip*.mat', fileName));
nDipMaps = size(dipMaps,1);

%% calculate fit for the region
fit = fit_source('filePath',filePath, 'expData', bData, ...
                'xy', [col(1), row(1)], 'dx', diff(col), 'dy', diff(row), ...
                'fitOrder', kwargs.fitOrder, 'nRuns', 2, ...
                'cropFactor', max([diff(col) diff(row)]), 'save', false,...
                'checkPlot', kwargs.checkPlot);

bData(row(1):row(2), col(1):col(2)) = fit.residuals;

expData.(dataName) = bData;

if kwargs.checkPlot
    QDM_figure(bData);
end

% calculate new Btotal
[by, bx] = MITBxByFromBz(bData, 1/expData.step);
expData.Bt = sqrt(by.^2+bx.^2+bData.^2);

if kwargs.save
    fit.datetime = datetime;
    expData.(sprintf('dipSub_%i', nDipMaps+1)) = fit;
    save(fullfile(path, [fileName, sprintf('-dip%i', nDipMaps+1), '.mat']), '-struct', 'expData');
end

end