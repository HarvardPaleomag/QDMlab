% This script takes an input Bz map, asks for a box, fits a source (D, Q,
% or O), and saves (1) a text file with the subtracted source parameters
% and (2) a new Bz map with the source subtracted. It prints a map of this
% new Bz map as well.  

%1=dipole; 2=quadrupole; 3=octapole
fitorder=1;


%Load the correct Bz file
[longfilename, pathname] = uigetfile('*.mat', 'Pick a magnetic field map file');
fullfilename=[pathname longfilename];
[filepath,name,ext]=fileparts(fullfilename);

load(fullfilename)

%Crop a single region for source fitting
disp(sprintf('Cropping area selection:  (+) saturate color scale, (-) desaturate color scale, \n(*) or (/) restore original color scale'));
figure;
imagesc(Bz);
title('Select area to crop')
axis equal, axis tight, axis xy
caxis([-1 1]*max(abs(caxis)));
caxis(caxis()/sqrt(10));
colormap(jet);
colorbar
flag=1;

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

residualmap=FitMoment(fitorder, fullfilename, [col(1,1) lin(1,1)], [col(2,1) lin(2,1)], 2, 0);
BzSub=Bz;

ic=1;
for i=col(1):col(2)
    jc=1;
    for j=lin(1):lin(2)
        BzSub(j,i)=-residualmap(jc,ic);
        jc=jc+1;
    end
    ic=ic+1;
end

figure
imagesc(BzSub)
caxis([-1 1]*max(abs(caxis)));
axis xy, axis equal, axis tight, axis off
hh=colorbar;
colormap(jet);

Bz=BzSub;
[by,bx]=MITBxByFromBz(Bz,1/step);
Bt=sqrt(by.^2+bx.^2+Bz.^2);

if size(str2num(name(end)),1)
    save([filepath '/' name(1:end-1) num2str(str2num(name(end))+1) '.mat'],'Bz','Bt','h','step','-mat');
else
    save([filepath '/' name '1' '.mat'],'Bz','Bt','h','step','-mat');
end





