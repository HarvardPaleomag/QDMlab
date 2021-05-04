function [] = B111ToBz_CommLine(INFILE,ferroorpara,pixel,h,quadsub,colorrange)
%This script takes in a B111 map file that must include a B111ferro or
%B111para matrix and an ledImg matrix
%for pixel size
%use 2.42 µm for QDM1 4x4 binning
%use 4.68 µm for the 10x objective with 4x4 binning
%use 2.34 µm for the 20x objective with 4x4 binning
%assumed 100 mm tube lens
%unlike the QDMfileviewer.m script, there's no cropping, interpolation, bkg subtraction,
%or upward continuation here.
%quadsub tells the script whether to apply a quadratic background fit
%leaving colorrange empty defaults to the max/min map values

UNIT='T';
SAVE=1;
ROTATE=0;
IFACTOR=1;
CAL=1;          %field calibration factor
alpha=0;        %rotation of the diamond lattice axes around z-axis
beta=0;         %rotation of the image axes around z-axis

[filepath,name,ext]=fileparts(INFILE);
load(INFILE);

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
        disp('Must specify "ferro" or "para" in second argument.');
    end
end

xstep=pixel;
h=h;

U=[ 0 sqrt(2/3) sqrt(1/3) ];%QDM2
%U=[ -sqrt(2/3) 0 -sqrt(1/3) ];%QDM1
disp(sprintf('Assumed U vector is:'));
disp(U);

Ra=[cosd(alpha) -sind(alpha)  0;
    sind(alpha) cos(alpha)    0
    0           0        1];


U=Ra*(U');
U_hat=U/norm(U);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ROTATE
    B=rot90(B,ROTATE);
end

FLIP=0;
Calib=CAL;

%scale and quadratic fit subtraction
if quadsub
    disp('Doing Quad BG subtraction...');
    Bback=QuadBGsub(B);
    s=Bback*Calib;
else
    s=B*Calib;
end

[bz]=QDMBzFromBu(s,1/xstep,U_hat);
[by,bx]=MITBxByFromBz(s,1/xstep);      %CHECK THIS IS RIGHT

bu=s;
bt=sqrt(bx.^2+by.^2+bz.^2);
Bz=bz;
Bu=bu;
Bt=bt;
step=xstep;
h=h;
newLED=ledImg;

corners=0;

if exist('corners')
    corners=corners;
else
    corners=[1,size(Bz,1);1,size(Bz,2)];
end

%show figures
figure
imagesc(Bz);
if isempty(colorrange)
    
    caxis([-1 1]*max(abs(caxis))*0.5);
else
    caxis(colorrange)
end
axis xy, axis equal, axis tight, axis off
hh=colorbar;
set(gca,'Fontsize',14);
title(hh,sprintf('       B_z (%s)',UNIT),'Fontsize',14);
colormap(jet);

saveas(gcf,[filepath '/BzCropped.png'])

if SAVE
    msg = sprintf('SAVING data to %s/Bz.mat', filepath');
    logMsg('info',msg,1,0);
    eval(['save ' filepath '/Bz.mat Bz Bt h step corners newLED -mat']);
end