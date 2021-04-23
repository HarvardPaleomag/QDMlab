%Load the correct "B111dataToPlot.mat" file
[longfilename, pathname] = uigetfile('*.mat', 'Pick a magnetic field map file');
fullfilename=[pathname longfilename];
[filepath,name,ext]=fileparts(fullfilename);

load(fullfilename)


%choose value for subtraction
disp(sprintf('Cropping area selection:  (+) saturate color scale, (-) desaturate color scale, \n(*) or (/) restore original color scale'));
figure;
imagesc(Bz);
title('Select constant to subtract')
axis equal, axis tight, axis xy
caxis([-1 1]*max(abs(caxis)));
colorbar
colormap(jet);
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
    flag=length(lin)~=1;
end
const=mean(mean(Bz(lin-2:lin+2,col-2:col+2)));
Bz=Bz-const;
figure;
imagesc(Bz);
title('Subtracted map')
axis equal, axis tight, axis xy
caxis([-1 1]*max(abs(caxis)));
colormap(jet);
colorbar
saveas(gcf,[filepath '/BzSub.png'])

disp(['Constant Bz value of ' num2str(const) ' subtracted from map']);

save([filepath '/BzSub.mat'],'Bz','Bt','h','step','corners','newLED');




