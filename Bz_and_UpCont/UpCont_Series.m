%function [] = QDMfileviewer()

UCarray=[10 20 30 40 50]*1e-6;%upward continuations in microns
%UCarray=[10 20 30 40 60 80 100]*1e-6;%upward continuations in microns
%UCarray=[60]*1e-6;%upward continuations in microns
%UCarray=[10,20,35,50,75,100,150]*1e-6;

UNIT='T';
SAVE=1;
ROTATE=0;
MEDFILTER=0;
CAL=1;          %field calibration factor
alpha=0;        %rotation of the diamond lattice axes around z-axis
beta=0;         %rotation of the image axes around z-axis
[longfilename, pathname] = uigetfile('*.mat', 'Pick a magnetic field map file');

fullfilename=[pathname longfilename];
parsedpath=strsplit(pathname,filesep);
originalfilename=extractBefore(longfilename,".mat");

load(fullfilename);

xstep=step;
%make even sized arrays
if mod(size(Bz,1),2)
    if mod(size(Bz,2),2)
        Bztruncated=Bz(2:end,2:end);
    else
        Bztruncated=Bz(2:end,1:end);
    end
else
    if mod(size(Bz,2),2)
        Bztruncated=Bz(1:end,2:end);
    else
        Bztruncated=Bz(1:end,1:end);
    end
end

updist=0;

[bz]=Bztruncated;

[by,bx]=MITBxByFromBz(bz,1/xstep);

%show figures
bt=sqrt(bx.^2+by.^2+bz.^2);
Bz=bz; %Conserve units
Bt=bt;
step=xstep;
h=h+updist;

figure
imagesc(Bz);
caxis([-1 1]*max(abs(caxis)));
axis xy, axis equal, axis tight, axis off
hh=colorbar;
set(gca,'Fontsize',14);
title(hh,sprintf('       B_z (%s)',UNIT),'Fontsize',14);
title('Original');
colormap(jet);

if SAVE
    save([pathname originalfilename '_uc' num2str(updist*1e6) '.mat'],'Bz', 'Bt', 'h', 'step', '-mat');
    disp(sprintf(['\nSaving Bz_uc' num2str(updist) '.mat...\n']))
end

for n=1:size(UCarray,2)
    updist=UCarray(n);
    
    [bz]=UpCont(Bztruncated,updist,1/(xstep));          %upward continue Bz map
    
    [by,bx]=MITBxByFromBz(bz,1/xstep);
    
    %show figures
    bt=sqrt(bx.^2+by.^2+bz.^2);
    Bz=bz; %Conserve units
    Bt=bt;
    step=xstep;
    h=h+updist;
    
    figure
    imagesc(Bz);
    caxis([-1 1]*max(abs(caxis)));
    axis xy, axis equal, axis tight, axis off
    hh=colorbar;
    set(gca,'Fontsize',14);
    title(hh,sprintf('       B_z (%s)',UNIT),'Fontsize',14);
    title(['Upward continuation ' num2str(updist*1e6) ' µm']);
    colormap(jet);
    
    if SAVE
        saveas(gcf,[pathname originalfilename '_uc' num2str(updist*1e6) '.png']);
        save([pathname originalfilename '_uc' num2str(updist*1e6) '.mat'],'Bz', 'Bt', 'h', 'step', '-mat');
        disp(sprintf(['\nSaving ' originalfilename '_uc' num2str(updist*1e6) '.mat...\n']))
    end
end
