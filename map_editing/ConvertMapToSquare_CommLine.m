function [] = ConvertMapToSquare_CommLine(INFILE)
% This script takes a Bz file and makes it square maintaining the center
% and preserving as much area as possible.
[filepath,name,ext]=fileparts(INFILE);
load(INFILE);
if size(Bz,1)==size(Bz,2)
    save([filepath '/' name '_SQ.mat'],'Bz','Bt','h','step','corners','newLED');
else
    binsize=round(size(newLED,1) / size(Bz,1));
    
    [row,colmax]=find(size(Bz)==max(size(Bz)));
    
    if colmax == 1
        sidecrop=(max(size(Bz))-min(size(Bz)));
        if ~mod(sidecrop,2)
            %cropping same amount each side
            low=sidecrop/2+1;
            high=size(Bz,1)-sidecrop/2;
            Bz=Bz(low:high,:);
            Bt=Bt(low:high,:);
            newLED=newLED(low*binsize:high*binsize,:);
            corners=[low,high;corners(2,1),corners(2,2)];
        else
            %cropping one pixel more on the lower side
            low=(sidecrop+1)/2+1;
            high=size(Bz,1)-(sidecrop-1)/2;
            Bz=Bz(low:high,:);
            Bt=Bt(low:high,:);
            newLED=newLED(low*binsize:high*binsize,:);
            corners=[low,high;corners(2,1),corners(2,2)];
        end
    else
        sidecrop=(max(size(Bz))-min(size(Bz)));
        if ~mod(sidecrop,2)
            %cropping same amount each side
            low=sidecrop/2+1;
            high=size(Bz,1)-sidecrop/2;
            Bz=Bz(:,low:high);
            Bt=Bt(:,low:high);
            newLED=newLED(:,low*binsize:high*binsize);
            corners=[corners(1,1),corners(1,2),low,high];
        else
            %cropping one pixel more on the lower side
            low=(sidecrop+1)/2+1;
            high=size(Bz,1)-(sidecrop-1)/2;
            Bz=Bz(:,low:high);
            Bt=Bt(:,low:high);
            newLED=newLED(:,low*binsize:high*binsize);
            corners=[corners(1,1),corners(1,2),low,high];
        end
    end
    figure
    imagesc(Bz)
    caxis([-1 1]*max(abs(caxis))*0.5);
    axis xy, axis equal, axis tight, axis off
    hh=colorbar;
    colormap(jet);
    saveas(gcf,[filepath '/' name '_SQ.png'])
    
    save([filepath '/' name '_SQ.mat'],'Bz','Bt','h','step','corners','newLED');
end





