function out = RMSInBox(kwargs)
% This script takes an input Bz map, asks for a box, crops to that box, and
% outputs Bz and Bt maps, along with the accessory parameters

arguments
    kwargs.dataFile = 'none'
    kwargs.binsize=4;
end

%Load the correct Bz file
% [longfilename, pathname] = uigetfile('*.mat', 'Pick a magnetic field map file');
% fullfilename=[pathname longfilename];
% [filepath,name,ext]=fileparts(fullfilename);
% 
% clear B111ferro;
% clear Bz;

dataFile = automatic_input_ui__(kwargs.dataFile, 'type', 'file', 'title', 'Pick a magnetic field map file');
expData = load(dataFile{:});
[~, dataName, ~] = is_B111(expData);

[lin, col] = pick_box2('expData', expData, 'title', 'Select area to crop');

%enforce that the cropped array has even dimensions
xrange=col(2)-col(1);
yrange=lin(2)-lin(1);
if mod(xrange,2)
    col(2,1)=col(2,1);
else
    col(2,1)=col(2,1)-1;
end
if mod(yrange,2)
    lin(2,1)=lin(2,1);
else
    lin(2,1)=lin(2,1)-1;
end


bData = expData.(dataName);
%cropping a B111 map set
bData=bData(lin(1):lin(2),col(1):col(2));
bData=bData-mean(mean(bData));

out = rms(rms(bData));

fig = figure('Units', 'normalized', ...
             'Position',[0.2 0.2 0.5 0.5], 'Name', sprintf('RMS: %.3e', out));
         
QDM_figure(bData, 'title', sprintf('cropped map RMS: %.3e', out), 'fig', fig)

fprintf('<>   The RMS of the selected region is: %.3e', out)




