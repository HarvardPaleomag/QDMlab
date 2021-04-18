function out = RMS_in_box(kwargs)
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
load(dataFile{:})

%Crop a single region for source fitting
if exist('B111ferro')
    disp(sprintf('Cropping area selection:  (+) saturate color scale, (-) desaturate color scale, \n(*) or (/) restore original color scale'));
    figure;
    imagesc(B111ferro);
    title('Select area to crop')
    axis equal, axis tight, axis xy
    caxis([-1 1]*max(abs(caxis)));
    caxis(caxis()/sqrt(10));
    colorbar
    colormap(jet);
    flag=1;
    
    binsize=1;
    if exist('ledImg')
        binsize=round(size(ledImg,1) / size(B111ferro,1));
    end
else
    figure;
    imagesc(Bz);
    title('Select area to crop')
    axis equal, axis tight, axis xy
    caxis([-1 1]*max(abs(caxis)));
    caxis(caxis()/sqrt(10));
    colorbar
    colormap(jet);
    flag=1;
    
    binsize=1;
    if exist('newLED')
        binsize=round(size(newLED,1) / size(Bz,1));
    end
end

lin=[];
col=[];
while flag
    [c,l,key]=ginput(1);
    switch key
        case 43
            caxis(caxis/sqrt(10));
        case 45
            caxis(caxis*sqrt(10));
        case {42,47}
            caxis auto, caxis([-1 1]*max(abs(caxis)));
        case 1
            lin=[lin; round(l)];
            col=[col; round(c)];
        case 27
            return
    end
    flag=length(lin)~=2;
end
lin=sort(lin,1);
col=sort(col,1);

if exist('B111ferro')
    if lin(1)<=0
        lin(1)=1;
    end
    if col(1)<=0
        col(1)=1;
    end
    if lin(2)>size(B111ferro,1)
        lin(2)=size(B111ferro,1);
    end
    if col(2)>size(B111ferro,2)
        col(2)=size(B111ferro,2);
    end
else
    if lin(1)<=0
        lin(1)=1;
    end
    if col(1)<=0
        col(1)=1;
    end
    if lin(2)>size(Bz,1)
        lin(2)=size(Bz,1);
    end
    if col(2)>size(Bz,2)
        col(2)=size(Bz,2);
    end
end

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

if exist('B111ferro')
    %cropping a B111 map set
    B111ferro=B111ferro(lin(1):lin(2),col(1):col(2));
    B111ferro=B111ferro-mean(mean(B111ferro));
    
    figure
    imagesc(B111ferro)
    caxis([-1 1]*max(abs(caxis))*0.5);
    axis xy, axis equal, axis tight, axis off
    hh=colorbar;
    colormap(jet);
    
    out = rms(rms(B111ferro));
else
    Bz=Bz(lin(1):lin(2),col(1):col(2));
    Bz=Bz-mean(mean(Bz));
    
    figure
    imagesc(Bz)
    caxis([-1 1]*max(abs(caxis))*0.5);
    axis xy, axis equal, axis tight, axis off
    hh=colorbar;
    colormap(jet);
    
    out = rms(rms(Bz));
end

fprintf('<>   The RMS of the selected region is: %.3e', out)




