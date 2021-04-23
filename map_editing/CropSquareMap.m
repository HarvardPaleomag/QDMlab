% This script takes an input B111 OR Bz map, asks for a box, makes a evenxeven square trying to maintain the same center point, and
% outputs cropped square Bz and Bt maps, along with the accessory parameters

binsize=1;
tempname='SQ';%leave as empty string if using default name. exclude extension

%Load the correct Bz file
[longfilename, pathname] = uigetfile('*.mat', 'Pick a magnetic field map file');
fullfilename=[pathname longfilename];
[filepath,name,ext]=fileparts(fullfilename);

clear B111ferro;
clear Bz;
clear ledImg;
clear newLED;
load(fullfilename)

%Crop a single region for source fitting
if exist('B111ferro')
    disp(sprintf('Cropping area selection:  (+) saturate color scale, (-) desaturate color scale, \n(*) or (/) restore original color scale'));
    figure;
    imagesc(B111ferro);
    title('Select box to crop')
    axis equal, axis tight, axis xy
    caxis([-1 1]*max(abs(caxis)));
    caxis(caxis()/sqrt(10));
    colorbar
    colormap(jet);
    flag=1;
    
    if exist('ledImg')
        binsize=round(size(ledImg,1) / size(B111ferro,1));
    end
else
    figure;
    imagesc(Bz);
    title('Select box to crop')
    axis equal, axis tight, axis xy
    caxis([-1 1]*max(abs(caxis)));
    caxis(caxis()/sqrt(10));
    colorbar
    colormap(jet);
    flag=1;
    
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

%make lin and col even
if mod(col(2)-col(1),2)
    col(2)=col(2)-1;
end
if mod(lin(2)-lin(1),2)
    lin(2)=lin(2)-1;
end

%make lin and col into a square
ww=col(2)-col(1);
hh=lin(2)-lin(1);
if ww > hh
    col(1)=(ww-hh)/2+col(1);
    col(2)=-(ww-hh)/2+col(2);
else
    if hh > ww
        lin(1)=(hh-ww)/2+lin(1);
        lin(2)=-(hh-ww)/2+lin(2);
    end
end

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
    B111para=B111para(lin(1):lin(2),col(1):col(2));
    if exist('ledImg')
        ledImg=ledImg(lin(1)*binsize:lin(2)*binsize,col(1)*binsize:col(2)*binsize);
    end
    corners=[lin(1),lin(2);col(1),col(2)];
    
    figure
    imagesc(B111ferro)
    caxis([-1 1]*max(abs(caxis))*0.5);
    axis xy, axis equal, axis tight, axis off
    hh=colorbar;
    colormap(jet);
    saveas(gcf,[filepath '/B111Cropped.png'])
    
    if exist('ledImg')
        if exist('corners')
            save([filepath '/B111_' num2str(rebin) 'x' num2str(rebin) 'Binned.mat'],'B111ferro','B111para','corners','ledImg');
        else
            save([filepath '/B111_' num2str(rebin) 'x' num2str(rebin) 'Binned.mat'],'B111ferro','B111para','ledImg');
        end
    else
        if exist('corners')
            save([filepath '/B111_' num2str(rebin) 'x' num2str(rebin) 'Binned.mat'],'B111ferro','B111para','corners');
        else
            save([filepath '/B111_' num2str(rebin) 'x' num2str(rebin) 'Binned.mat'],'B111ferro','B111para');
        end
    end
else
    Bz=Bz(lin(1):lin(2),col(1):col(2));
    Bt=Bt(lin(1):lin(2),col(1):col(2));
    if exist('newLED')
        newLED=newLED(lin(1)*binsize:lin(2)*binsize,col(1)*binsize:col(2)*binsize);
    end
    corners=[lin(1),lin(2);col(1),col(2)];
    
    figure
    imagesc(Bz)
    caxis([-1 1]*max(abs(caxis))*0.5);
    axis xy, axis equal, axis tight, axis off
    hh=colorbar;
    colormap(jet);
    if isempty(tempname)
        saveas(gcf,[filepath '/' name '_SQ.png'])
    else
        saveas(gcf,[filepath '/' tempname '.png'])
    end
    
    if isempty(tempname)
        if exist('newLED')
            if exist('corners')
                save([filepath '/' name '_SQ.mat'],'Bz','Bt','h','step','corners','newLED');
            else
                save([filepath '/' name '_SQ.mat'],'Bz','Bt','h','step','newLED');
            end
        else
            if exist('corners')
                save([filepath '/' name '_SQ.mat'],'Bz','Bt','h','step','corners');
            else
                save([filepath '/' name '_SQ.mat'],'Bz','Bt','h','step');
            end
        end
    else
        if exist('newLED')
            if exist('corners')
                save([filepath '/' tempname '.mat'],'Bz','Bt','h','step','corners','newLED');
            else
                save([filepath '/' tempname '.mat'],'Bz','Bt','h','step','newLED');
            end
        else
            if exist('corners')
                save([filepath '/' tempname '.mat'],'Bz','Bt','h','step','corners');
            else
                save([filepath '/' tempname '.mat'],'Bz','Bt','h','step');
            end
        end
    end
end


