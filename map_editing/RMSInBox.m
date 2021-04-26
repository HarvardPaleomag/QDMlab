function out = RMSInBox(kwargs)
% This script takes an input Bz map, asks for a box, crops to that box, and
% outputs Bz and Bt maps, along with the accessory parameters

arguments
    kwargs.dataFile = 'none'
    kwargs.binsize=4;
    kwargs.saveCropMap=false;
end

expData = crop_map('dataFile', kwargs.dataFile, 'cropFigure', false, ...
                'save', kwargs.saveCropMap);
[~, dataName, ~] = is_B111(expData);
bData = expData.(dataName);
bData = bData-mean(mean(bData));

out = rms(rms(bData));

fig = figure('Units', 'normalized', ...
             'Position',[0.2 0.2 0.5 0.5], 'Name', sprintf('RMS: %.3e', out));
         
QDM_figure(bData, 'title', sprintf('cropped map RMS: %.3e', out), 'fig', fig);

fprintf('<>   The RMS of the selected region is: %.3e', out);




