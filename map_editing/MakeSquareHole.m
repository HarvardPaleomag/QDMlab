[longfilename, pathname] = uigetfile('*.mat', 'Pick a magnetic field map file');
fullfilename=[pathname longfilename];
parsedpath=strsplit(pathname,filesep);
originalfilename=extractBefore(longfilename,".mat");

ferroorpara=input('[ferro] or [para] map? [ferro]: ','s');
if isempty(ferroorpara)
    ferroorpara = 'ferro';
end

if ferroorpara(1:4) =='para'
    ferroorpara ='para ';
end

if ferroorpara == 'ferro'
    clear B111ferro;
else
    clear B111para;
end
clear Bz;

load(fullfilename)

%Select a background value
if ferroorpara == 'ferro'
    if exist('B111ferro')
        disp(sprintf('Cropping area selection:  (+) saturate color scale, (-) desaturate color scale, \n(*) or (/) restore original color scale'));
        figure;
        imagesc(B111ferro);
        title('Select background value')
        axis equal, axis tight, axis xy
        caxis([-1 1]*max(abs(caxis)));
        caxis(caxis()/sqrt(10));
        colorbar
        colormap(jet);
        flag=1;
        
        flag=1;
        while flag
            [bkgx,bkgy,key] = ginput(1);
            switch key
                case 43
                    caxis(caxis/sqrt(10));
                case 45
                    caxis(caxis*sqrt(10));
                case 1
                    flag=0;
                    bkgx=round(bkgx);
                    bkgy=round(bkgy);
            end
        end
        
        limx=max(1,bkgx-2);
        limX=min(size(B111ferro,2),bkgx+2);
        limy=max(1,bkgy-2);
        limY=min(size(B111ferro,1),bkgy+2);
        bkg=mean2(B111ferro(limy:limY,limx:limX));
    else
        figure;
        imagesc(Bz);
        title('Select background value')
        axis equal, axis tight, axis xy
        caxis([-1 1]*max(abs(caxis)));
        caxis(caxis()/sqrt(10));
        colorbar
        colormap(jet);
        flag=1;
        
        flag=1;
        while flag
            [bkgx,bkgy,key] = ginput(1);
            switch key
                case 43
                    caxis(caxis/sqrt(10));
                case 45
                    caxis(caxis*sqrt(10));
                case 1
                    flag=0;
                    bkgx=round(bkgx);
                    bkgy=round(bkgy);
            end
        end
        
        limx=max(1,bkgx-2);
        limX=min(size(Bz,2),bkgx+2);
        limy=max(1,bkgy-2);
        limY=min(size(Bz,1),bkgy+2);
        bkg=mean2(Bz(limy:limY,limx:limX));
    end
    
    %choose the diagonal corners of the area to be removed
    if exist('B111ferro')
        disp(sprintf('Cropping area selection:  (+) saturate color scale, (-) desaturate color scale, \n(*) or (/) restore original color scale'));
        figure;
        imagesc(B111ferro);
        title('Select area to remove')
        axis equal, axis tight, axis xy
        caxis([-1 1]*max(abs(caxis)));
        caxis(caxis()/sqrt(10));
        colorbar
        colormap(jet);
        flag=1;
        
        xx=[];
        yy=[];
        flag=1;
        while flag
            [x,y,key]=ginput(1);
            switch key
                case 43
                    caxis(caxis/sqrt(10));
                case 45
                    caxis(caxis*sqrt(10));
                case 1
                    xx=[xx; round(x)];
                    yy=[yy; round(y)];
            end
            flag=length(xx)~=2;
        end
        
        Nreplacements=0;
        for i = min(xx):max(xx)
            for j = min(yy):max(yy)
                B111ferro(j,i) = bkg;
                Nreplacements = Nreplacements + 1;
            end
        end
        
        figure
        imagesc(B111ferro)
        caxis([-1 1]*max(abs(caxis)));
        axis xy, axis equal, axis tight, axis off
        hh=colorbar;
        colormap(jet);
        disp([num2str(Nreplacements) ' entries in B111ferro replaced by ' num2str(bkg)])
        
        saveas(gcf,[pathname originalfilename '_Hole.png'])
        save([pathname originalfilename '_Hole.mat'],'B111ferro', 'B111para', 'ledImg', 'negDiff', 'posDiff', '-mat')
        disp(sprintf(['\nSaving ' pathname originalfilename '_Hole.mat...\n']))
    else
        figure;
        imagesc(Bz);
        title('Select area to remove')
        axis equal, axis tight, axis xy
        caxis([-1 1]*max(abs(caxis)));
        caxis(caxis()/sqrt(10));
        colorbar
        colormap(jet);
        flag=1;
        
        xx=[];
        yy=[];
        flag=1;
        while flag
            [x,y,key]=ginput(1);
            switch key
                case 43
                    caxis(caxis/sqrt(10));
                case 45
                    caxis(caxis*sqrt(10));
                case 1
                    xx=[xx; round(x)];
                    yy=[yy; round(y)];
            end
            flag=length(xx)~=2;
        end
        
        Nreplacements=0;
        for i = min(xx):max(xx)
            for j = min(yy):max(yy)
                Bz(j,i) = bkg;
                Nreplacements = Nreplacements + 1;
            end
        end
        
        figure
        imagesc(Bz)
        caxis([-1 1]*max(abs(caxis)));
        axis xy, axis equal, axis tight, axis off
        hh=colorbar;
        colormap(jet);
        disp([num2str(Nreplacements) ' entries in Bz replaced by ' num2str(bkg)])
        
        saveas(gcf,[pathname originalfilename '_Hole.png'])
        save([pathname originalfilename '_Hole.mat'], 'Bt', 'Bz', 'corners', 'h', 'newLED', 'step', '-mat')
        disp(sprintf(['\nSaving ' pathname originalfilename '_Hole.mat...\n']))
    end
    
else
    if exist('B111para')
        
        plottablepara=B111para-mean(mean(B111para));
        
        disp(sprintf('Cropping area selection:  (+) saturate color scale, (-) desaturate color scale, \n(*) or (/) restore original color scale'));
        figure;
        imagesc(plottablepara);
        title('Select background value')
        axis equal, axis tight, axis xy
        caxis([-1 1]*max(abs(caxis)));
        caxis(caxis()/sqrt(10));
        colorbar
        colormap(jet);
        flag=1;
        
        flag=1;
        while flag
            [bkgx,bkgy,key] = ginput(1);
            switch key
                case 43
                    caxis(caxis/sqrt(10));
                case 45
                    caxis(caxis*sqrt(10));
                case 1
                    flag=0;
                    bkgx=round(bkgx);
                    bkgy=round(bkgy);
            end
        end
        
        limx=max(1,bkgx-2);
        limX=min(size(B111para,2),bkgx+2);
        limy=max(1,bkgy-2);
        limY=min(size(B111para,1),bkgy+2);
        bkg=mean2(B111para(limy:limY,limx:limX));
    else
        figure;
        imagesc(Bz);
        title('Select background value')
        axis equal, axis tight, axis xy
        caxis([-1 1]*max(abs(caxis)));
        caxis(caxis()/sqrt(10));
        colorbar
        colormap(jet);
        flag=1;
        
        flag=1;
        while flag
            [bkgx,bkgy,key] = ginput(1);
            switch key
                case 43
                    caxis(caxis/sqrt(10));
                case 45
                    caxis(caxis*sqrt(10));
                case 1
                    flag=0;
                    bkgx=round(bkgx);
                    bkgy=round(bkgy);
            end
        end
        
        limx=max(1,bkgx-2);
        limX=min(size(Bz,2),bkgx+2);
        limy=max(1,bkgy-2);
        limY=min(size(Bz,1),bkgy+2);
        bkg=mean2(Bz(limy:limY,limx:limX));
    end
    
    %choose the diagonal corners of the area to be removed
    if exist('B111para')
        disp(sprintf('Cropping area selection:  (+) saturate color scale, (-) desaturate color scale, \n(*) or (/) restore original color scale'));
        figure;
        imagesc(plottablepara);
        title('Select area to remove')
        axis equal, axis tight, axis xy
        caxis(mean(mean(B111para))+[-1 1]*max(abs(caxis)));
        caxis(caxis()/sqrt(10));
        colorbar
        colormap(jet);
        flag=1;
        
        xx=[];
        yy=[];
        flag=1;
        while flag
            [x,y,key]=ginput(1);
            switch key
                case 43
                    caxis(caxis/sqrt(10));
                case 45
                    caxis(caxis*sqrt(10));
                case 1
                    xx=[xx; round(x)];
                    yy=[yy; round(y)];
            end
            flag=length(xx)~=2;
        end
        
        Nreplacements=0;
        for i = min(xx):max(xx)
            for j = min(yy):max(yy)
                B111para(j,i) = bkg;
                Nreplacements = Nreplacements + 1;
            end
        end
        
        figure
        imagesc(B111para)
        caxis([-1 1]*max(abs(caxis)));
        axis xy, axis equal, axis tight, axis off
        hh=colorbar;
        colormap(jet);
        disp([num2str(Nreplacements) ' entries in B111para replaced by ' num2str(bkg)])
        
        saveas(gcf,[pathname originalfilename '_Hole.png'])
        save([pathname originalfilename '_Hole.mat'],'B111ferro', 'B111para', 'ledImg', 'negDiff', 'posDiff', '-mat')
        disp(sprintf(['\nSaving ' pathname originalfilename '_Hole.mat...\n']))
    else
        figure;
        imagesc(Bz);
        title('Select area to remove')
        axis equal, axis tight, axis xy
        caxis([-1 1]*max(abs(caxis)));
        caxis(caxis()/sqrt(10));
        colorbar
        colormap(jet);
        flag=1;
        
        xx=[];
        yy=[];
        flag=1;
        while flag
            [x,y,key]=ginput(1);
            switch key
                case 43
                    caxis(caxis/sqrt(10));
                case 45
                    caxis(caxis*sqrt(10));
                case 1
                    xx=[xx; round(x)];
                    yy=[yy; round(y)];
            end
            flag=length(xx)~=2;
        end
        
        Nreplacements=0;
        for i = min(xx):max(xx)
            for j = min(yy):max(yy)
                Bz(j,i) = bkg;
                Nreplacements = Nreplacements + 1;
            end
        end
        
        figure
        imagesc(Bz)
        caxis([-1 1]*max(abs(caxis)));
        axis xy, axis equal, axis tight, axis off
        hh=colorbar;
        colormap(jet);
        disp([num2str(Nreplacements) ' entries in Bz replaced by ' num2str(bkg)])
        
        saveas(gcf,[pathname originalfilename '_Hole.png'])
        save([pathname originalfilename '_Hole.mat'], 'Bt', 'Bz', 'corners', 'h', 'newLED', 'step', '-mat')
        disp(sprintf(['\nSaving ' pathname originalfilename '_Hole.mat...\n']))
    end
end


