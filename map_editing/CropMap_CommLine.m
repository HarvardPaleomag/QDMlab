function [] = CropMap_CommLine(INFILE,outname,center,sizes,colorrange)
% This script takes as input:
% (1) B111 or Bz map
% (2) center point
% (3) x and y crop dimensions in units of pixels
% For B111 maps, outputs are cropped B111ferro and B111para and LED maps
% For Bz maps outputs are Bz and Bt maps, along with the accessory parameters and LED
% map

xbounds=[center(1)-round(sizes(1)/2),center(1)+round(sizes(1)/2)];%these are columns
ybounds=[center(2)-round(sizes(2)/2),center(2)+round(sizes(2)/2)];%these are rows

%enforce that the cropped array has even dimensions
xrange=xbounds(2)-xbounds(1);
yrange=ybounds(2)-ybounds(1);
if mod(xrange,2)
    xbounds(2)=xbounds(2);
else
    xbounds(2)=xbounds(2)-1;
end
if mod(yrange,2)
    ybounds(2)=ybounds(2);
else
    ybounds(2)=ybounds(2)-1;
end

%Load the correct B111 or Bz file
[filepath,name,ext]=fileparts(INFILE);
disp(filepath)
step=0;
load(INFILE);


if exist('B111ferro')
    %cropping a B111 map set
    binsize=round(size(ledImg,1) / size(B111ferro,1));
    B111ferro=B111ferro(ybounds(1):ybounds(2),xbounds(1):xbounds(2));
    B111para=B111para(ybounds(1):ybounds(2),xbounds(1):xbounds(2));
    ledImg=ledImg(ybounds(1)*binsize:ybounds(2)*binsize,xbounds(1)*binsize:xbounds(2)*binsize);
    corners=[ybounds(1),ybounds(2);xbounds(1),xbounds(2)];
    
    figure
    imagesc(B111ferro)
    if isempty(colorrange)
        caxis([-1 1]*max(abs(caxis))*0.5);
    else
        caxis(colorrange)
    end
    axis xy, axis equal, axis tight, axis off
    hh=colorbar;
    colormap(jet);
    saveas(gcf,[filepath '/' outname '.png'])
    
    if exist('ledImg','var')
        save([filepath '/' outname '.mat'],'B111ferro','B111para','corners','ledImg');
    else
        save([filepath '/' outname '.mat'],'B111ferro','B111para');
    end
else
    Bz=Bz(ybounds(1):ybounds(2),xbounds(1):xbounds(2));
    if exist('newLED','var')
        binsize=round(size(newLED,1) / size(Bz,1));
        newLED=newLED(ybounds(1)*binsize:ybounds(2)*binsize,xbounds(1)*binsize:xbounds(2)*binsize);
        Bt=Bt(ybounds(1):ybounds(2),xbounds(1):xbounds(2));
        corners=[ybounds(1),ybounds(2);xbounds(1),xbounds(2)];
    end
    
    figure
    imagesc(Bz)
    if isempty(colorrange)
        caxis([-1 1]*max(abs(caxis))*0.5);
    else
        caxis(colorrange)
    end
    axis xy, axis equal, axis tight, axis off
    hh=colorbar;
    colormap(jet);
    saveas(gcf,[filepath '/' outname '.png'])
    
    if exist('newLED','var')
        save([filepath '/' outname '.mat'],'Bz','Bt','h','step','corners','newLED');
    else
        disp('hey');
        save([filepath '/' outname '.mat'],'Bz','h','step');
    end
end





