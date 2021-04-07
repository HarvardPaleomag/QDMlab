%Load the correct "B111dataToPlot.mat" file
[longfilename, pathname] = uigetfile('*.mat', 'Pick a magnetic field map file');
fullfilename=[pathname longfilename];
[filepath,name,ext]=fileparts(fullfilename);

load(fullfilename);

%choose ferro or para map
ferroorpara=input('[ferro] or [para] map? [ferro]: ','s');
if isempty(ferroorpara)
    ferroorpara = 'ferro';
end

if strcmp(ferroorpara,'ferro')
    %choose value for subtraction
    disp(sprintf('Cropping area selection:  (+) saturate color scale, (-) desaturate color scale, \n(*) or (/) restore original color scale'));
    figure;
    imagesc(B111ferro);
    title('Select constant to subtract')
    axis equal, axis tight, axis xy
    %caxis([-1 1]*max(abs(caxis)));
    caxis([-.1 .1]);
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
        flag=length(lin)~=1;
    end
    const=B111ferro(lin,col);
    B111ferro=B111ferro-const;
    figure;
    imagesc(B111ferro);
    title('Subtracted map')
    axis equal, axis tight, axis xy
    caxis([-1 1]*max(abs(caxis)));
    colorbar
    
    save([filepath '/B111dataToPlotSub.mat'],'negDiff','posDiff', 'B111ferro', 'B111para', 'ledImg');
else
        %choose value for subtraction
    disp(sprintf('Cropping area selection:  (+) saturate color scale, (-) desaturate color scale, \n(*) or (/) restore original color scale'));
    figure;
    imagesc(B111para);
    title('Select constant to subtract')
    axis equal, axis tight, axis xy
    caxis([min(abs(caxis)) max(abs(caxis))]);
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
        flag=length(lin)~=1;
    end
    const=B111para(lin,col);
    B111para=B111para-const;
    figure;
    imagesc(B111para);
    title('Subtracted map')
    axis equal, axis tight, axis xy
    caxis([-1 1]*max(abs(caxis)));
    colorbar
    
    save([filepath '/B111dataToPlotSub.mat'],'negDiff','posDiff', 'B111ferro', 'B111para', 'ledImg');
end



