function out = RMSInBox(kwargs)
%[out] = RMSInBox('binsize', 'filePath', 'saveCropMap')
% This script takes an input Bz map, asks for a box, crops to that box, and
% outputs Bz and Bt maps, along with the accessory parameters

arguments
    kwargs.filePath = 'none'
    kwargs.binsize = 4;
    kwargs.saveCropMap = false;
end

expData = crop_map('filePath', kwargs.filePath, 'checkPlot', false, ...
    'save', kwargs.saveCropMap);

[~, dataName, ~] = is_B111(expData);
bData = expData.(dataName);
bData = bData - mean(mean(bData));

msg = sprintf('calculating RMS in the cropped region');
logMsg('info', msg, 1, 0);
out = rms(rms(bData));

fig = figure('Units', 'normalized', ...
    'Position', [0.2, 0.2, 0.5, 0.5], 'Name', sprintf('RMS: %.3e', out));

QDM_figure(bData, 'title', sprintf('cropped map RMS: %.3e', out), 'fig', fig);

msg = sprintf('The RMS of the selected region is: %.3e', out);
logMsg('Result', msg, 1, 1);
