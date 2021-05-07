% QDMdataprocessing.m -- MATLAB script for performing cropping, background compensation,
%   interpolation, calibration, and upward continuation of QDM data.
%   Processed data are saved to a .mat file and images are saved to .png
%   files.
%
%  Reference: R.R. Fu, E.A. Lima, M.W.R. Volk, R. Trubko (2020) "High-sensitivity
%  moment magnetometry with the quantum diamond microscope,"  Geochemistry,
%  Geophysics, Geosystems.
%
% -----------------------------------------------------------------------------------
%                Matlab code by Eduardo A. Lima and Roger R. Fu
%      Copyright (C) 2017-2020 MIT Paleomagnetism Lab, Harvard Paleomagnetics Lab
% -----------------------------------------------------------------------------------

show_references()

UNIT='T';
SAVE=1;         %1 = save processed data to disk
MEDFILTER=0;    %1 = apply median filter to reduce noise
CAL=1;          %field calibration factor
alpha=0;        %rotation of the diamond lattice axes around z-axis
beta=0;         %rotation of the image axes around z-axis

set(0,'DefaultFigureColormap',jet)      % Comment this if you prefer parula color map

%select and load data file
[longfilename, pathname] = uigetfile('*.mat', 'Pick a magnetic field map file');
fullfilename=[pathname longfilename];
parsedpath=strsplit(pathname,filesep);
savefilename=parsedpath(end-1);
savefilename=char(savefilename);

load(fullfilename);

%rotate image and adjust field units
ferroorpara=input('[ferro] or [para] map? [ferro]: ','s');
if isempty(ferroorpara)
    ferroorpara = 'ferro';
end

if strcmp(ferroorpara,'ferro')
    if beta
        B=imrotate(B111ferro,beta,'bilinear') * 1e-4; % Field in T
    else
        B=imrotate(B111ferro,0) * 1e-4; % Field in T
    end
else
    if strcmp(ferroorpara,'para')
        if beta
            B=imrotate(B111para,beta,'bilinear') * 1e-4 ; % Field in T
        else
            B=imrotate(B111para,0) * 1e-4; % Field in T
        end
    else
        disp('Must specify "ferro" or "para". ');
        return
    end
end

%ask user for the pixel size
%use 2.42 µm for QDM1 4x4 binning
%use 4.68 µm for the 10x objective with 4x4 binning
%use 2.34 µm for the 20x objective with 4x4 binning
%assumed 100 mm tube lens

xstep=input('Pixel size [4.70 µm]: ')*1e-6;
if isempty(xstep)
    xstep=4.70e-6;
end

%ask user for the NV layer-sample distance
h=input('NV-sample distance [5 µm]: ');
if isempty(h)
    h=5e-6;
end

%Unit vector specifying B111 direction
U=[ 0 sqrt(2/3) sqrt(1/3) ];%QDM2
%U=[ -sqrt(2/3) 0 -sqrt(1/3) ];%QDM1
disp(sprintf('Assumed U vector is:'));
disp(U);

%Apply rotation to unit vector if necessary
Ra=[cosd(alpha) -sind(alpha)  0;
    sind(alpha) cos(alpha)    0
    0           0        1];


U=Ra*(U');
U_hat=U/norm(U);

LEDcropfactor=size(ledImg,1) / size(B111ferro,1);

%tools for selecting a region of interest and cropping the field map
%and optical image accordingly
lin=[];
col=[];
disp(sprintf('Cropping area selection:  (+) saturate color scale, (-) desaturate color scale, \n(*) or (/) restore original color scale, (S)kip'));
figure;
imagesc(B);
title('Select area to crop')
axis equal, axis tight, axis xy
caxis([-1 1]*max(abs(caxis))*0.1);
colorbar
flag=1;
while flag
    [c,l,key]=ginput(1);
    switch key
        case 43           %saturate color scale (applied twice saturates by a factor of 10)
            caxis(caxis/sqrt(10));
        case 45           %desaturate color scale
            caxis(caxis*sqrt(10));
        case {42,47}      %restore normal color scale
            caxis auto, caxis([-1 1]*max(abs(caxis)));
        case 1
            lin=[lin; round(l)];    %select one corner of the cropping region
            col=[col; round(c)];
        case 27                     %quit code if user hits Esc
            return
    end
    if key=='S' || key=='s'          %skip cropping
        lin=[1; size(B,1)];
        col=[1; size(B,2)];
        flag=0;
    else
        flag=length(lin)~=2;
    end
end

lin=sort(lin,1);
col=sort(col,1);

%make sure cropping boundaries do not extend past image boundaries
for nn=1:2
    if lin(nn)<1
        lin(nn)=1;
    end
    if lin(nn)>size(B,1)
        lin(nn)=size(B,1);
    end
    if col(nn)<1
        col(nn)=1;
    end
    if col(nn)>size(B,2)
        col(nn)=size(B,2);
    end
end
%plot magenta rectangle showing the cropping area
hold on
plot([col(1) col(2) col(2) col(1) col(1)],[lin(1) lin(1) lin(2) lin(2) lin(1)],'m--');
hold off


corners=[lin,col];

B=B(max([1 lin(1)]):min([lin(2) size(B,1)]) , max([1 col(1)]):min([col(2) size(B,2)]) );

%crop the LED image in the same way
croppoint1=round(1+(corners(1,1)-1)*LEDcropfactor);
croppoint2=round(corners(2,1)*LEDcropfactor);
croppoint3=round(1+(corners(1,2)-1)*LEDcropfactor);
croppoint4=round(corners(2,2)*LEDcropfactor);

newLED=ledImg(croppoint1:croppoint2,croppoint3:croppoint4);

%discard one row or one column if necessary so as to make B have an even number of points
%in each coordinate
if mod(size(B,2),2)
    B=B(:,2:end);
end
if mod(size(B,1),2)
    B=B(2:end,:);
end
close(1);

disp('Background removal  - (R)evert (C)onstant (D)one: ');
figure(1);
flag=1;
Bback=B; %Bback is the background subtracted B111 map
ca=[-1 1]*max(abs(Bback(:)));  %ca is the default color scale
while flag
    imagesc(Bback);
    title('Background removal');
    axis equal, axis tight, axis xy
    caxis(ca);
    waitforbuttonpress;
    key=upper(get(gcf,'currentcharacter'));
    
    switch key
        case 43             %saturate color scale
            ca=ca/sqrt(10);
        case 45             %desaturate color scale
            ca=ca*sqrt(10);
        case {42,47}        %restore color scale
            ca=[-1 1]*max(abs(Bback(:)));
        case 27
            return
        case 'C'
            [c,l,kk]=ginput(1); % select center of 5 x 5 reference background region
            if kk==27
                return
            end
            lin=min([max([round(l),1]) size(B,1)]);
            col=min([max([round(c),1]) size(B,2)]);
        case 'R'            %undo any background compensation and revert to original data
            Bback=B;
            disp('Reverting to original data...');
        case 'D'            %set flag to exit backgroung compensation tool
            flag=0;
            
    end
    
    switch key
        case 'C'
            limx=max(1,col-2);      %make sure selected background region does not extend beyond boundaries
            limX=min(size(B,2),col+2);
            limy=max(1,lin-2);
            limY=min(size(B,1),lin+2);
            bkg=mean2(Bback(limy:limY,limx:limX));
            Bback=Bback-bkg;        %substract mean background value
            disp('Subtracting average field value...')
            disp(['Subtracted value is:' num2str(bkg)])
    end
end

Calib=CAL;

%perform quadratic background subtraction
quadbkg=input('Perfom quadratic fitting and subtraction [N]:','s');

if quadbkg=='y' | quadbkg=='Y'
    disp('Doing Quad BG subtraction...');
    Bback=QuadBGsub(Bback);
end

s=Bback*Calib;

%ask user for the interpolation factor to be applied to the data (1 = no
%interpolation)
IntFact=input('Interpolation factor [1]): ','s');
if isempty(IntFact)
    IFACTOR=1;
else
    IFACTOR=str2num(IntFact);
end

%select upward continuation distance (0 = no upward continuation)
up=input('Upward continuation distance [0]): ','s');
if isempty(up)
    updist=0e-6;
else
    updist=str2double(up)*1e-6;
end

sspreup=s;

if updist==0
    ss=sspreup;     %ss is the B111 map after cropping, bkg sub, and upcont
else
    ss=UpCont(sspreup,updist,1/xstep);          %upward continue B111 map
end

[bz]=QDMBzFromBu(ss,1/xstep,U_hat);             %compute Bz from B111
[by,bx]=MITBxByFromBz(ss,1/xstep);              %compute Bx and By from Bz 
bu=ss;

%calculate interpolated field maps
xm=0;
xM=size(B,2)-1;
ym=0;
yM=size(B,1)-1;
xx=xm:(xM-xm)/(size(B,2)-1):xM;
yy=ym:(yM-ym)/(size(B,1)-1):yM;
xxi=xm:(xM-xm)/(size(B,2)*IFACTOR-1):xM;
xyi=ym:(yM-ym)/(size(B,1)*IFACTOR-1):yM;
[XX,YY]=meshgrid(xx,yy);
[XXi,YYi]=meshgrid(xxi,xyi);
bzi=interp2(XX,YY,bz,XXi,YYi,'linear');
bui=interp2(XX,YY,bu,XXi,YYi,'linear');
byi=interp2(XX,YY,by,XXi,YYi,'linear');
bxi=interp2(XX,YY,bx,XXi,YYi,'linear');

%apply 4 x 4 median filter to Bz and Bu, if selected
if MEDFILTER
    bzi=medfilt2(bzi,[4 4]);
    bui=medfilt2(bui,[4 4]);
end


%show figures
figure
imagesc(bui);       % plot processed B111 map
caxis([-1 1]*max(abs(caxis)));
axis xy, axis equal, axis tight, axis off
hh=colorbar;
set(gca,'Fontsize',14);
title(hh,sprintf('       B_u (%s)',UNIT),'Fontsize',14);
colormap(jet);

bti=sqrt(bxi.^2+byi.^2+bzi.^2); %compute the total field
Bz=double(bzi);
Bu=double(bui);
Bt=double(bti);
step=xstep/IFACTOR;     %adjust step size to reflect interpolation
h=h+updist;             %adjust height to reflect upward continuation

figure
imagesc(Bz);        % plot calculated Bz map
caxis([-1 1]*max(abs(caxis)));
axis xy, axis equal, axis tight, axis off
hh=colorbar;
set(gca,'Fontsize',14);
title(hh,sprintf('       B_z (%s)',UNIT),'Fontsize',14);
colormap(jet);

%save images of (processed) B111 and Bz to disk
if length(ferroorpara)==4
    ferroorpara = [ferroorpara ' '];
end
if ferroorpara == 'ferro'
    saveas(gcf,[pathname 'Bz_uc' num2str(updist*1e6) '.png'])
end
if ferroorpara == 'para '
    saveas(gcf,[pathname 'ParaBz_uc' num2str(updist*1e6) '.png'])
end
caxis([-1 1]*2e-6);
if length(ferroorpara)==4
    ferroorpara = [ferroorpara ' '];
end
if ferroorpara == 'ferro'
    saveas(gcf,[pathname 'Bz_uc' num2str(updist*1e6) '_sat.png'])
end
if ferroorpara == 'para '
    saveas(gcf,[pathname 'ParaBz_uc' num2str(updist*1e6) '_sat.png'])
end

%save processed/calculated field maps, processed optical image, and ancillary info to disk
if SAVE
    if length(ferroorpara)==4
        ferroorpara = [ferroorpara ' '];
    end
    if ferroorpara == 'ferro'
        save([pathname 'Bz_uc' num2str(updist*1e6) '.mat'], 'Bz', 'Bt', 'h', 'step', 'corners', 'newLED', '-mat');
        disp(sprintf(['\nSaving ' savefilename '_Bz_uc' num2str(updist*1e6) '.mat...\n']))
    end
    if ferroorpara == 'para '
        save([pathname 'ParaBz_uc' num2str(updist*1e6) '.mat'], 'Bz', 'Bt', 'h', 'step', 'corners', 'newLED', '-mat');
        disp(sprintf(['\nSaving ' savefilename '_ParaBz_uc' num2str(updist*1e6) '.mat...\n']))
    end
end