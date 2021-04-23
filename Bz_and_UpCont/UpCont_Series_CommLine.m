function [] = UpCont_Series_CommLine(INFILE,UCarrayMICRONS)
%Input file must be a .mat file with a "Bz" variable that is a map with a
%"step" variable in METERS

if size(UCarrayMICRONS,2)~=0
    UCarray=UCarrayMICRONS*1e-6;%upward continuations in microns
end

UNIT='T';
SAVE=1;
ROTATE=0;
MEDFILTER=0;
CAL=1;          %field calibration factor
alpha=0;        %rotation of the diamond lattice axes around z-axis
beta=0;         %rotation of the image axes around z-axis
[filepath,name,ext] = fileparts(INFILE);

load(INFILE);

% if a B111dataToPlot.mat file is passed, the step, h variables is not 
% saved -> need to define it. 
if exist('step','var') == 0
    step=4.68e-6; % pixel size in (m)
    h = 5e-6;     % NV layer thickness (m) ?
end
    
% check filename for B111 data
% renamed the Bz variable to B_data
if name == 'B111dataToPlot'
    B_data = B111ferro;
else
    B_data = Bz;
end

xstep=step; %todo why rename the variable?

%make even sized arrays
if mod(size(B_data,1),2)
    if mod(size(B_data,2),2)
        B_truncated=B_data(2:end,2:end);
    else
        B_truncated=B_data(2:end,1:end);
    end
else
    if mod(size(B_data,2),2)
        B_truncated=B_data(1:end,2:end);
    else
        B_truncated=B_data(1:end,1:end);
    end
end

if size(UCarrayMICRONS,2)~=0
    for n=1:size(UCarray,2)
        updist=UCarray(n);

        if updist==0
            [b_uc]=B_truncated;
        else
            [b_uc]=UpCont(B_truncated,updist,1/(xstep)); %upward continue Bz map
        end
        
        [by,bx]=MITBxByFromBz(b_uc,1/xstep);
        
        %show figures
        bt=sqrt(bx.^2+by.^2+b_uc.^2);
        B_data=b_uc; %Conserve units
        Bt=bt;
        step=xstep; % rename variable again?
        h=h+updist;
        
        figure
        imagesc(B_data);
        caxis([-1 1]*max(abs(caxis)));
        axis xy, axis equal, axis tight, axis off
        hh=colorbar;
        set(gca,'Fontsize',14);
        title(hh,sprintf('       B_z (%s)',UNIT),'Fontsize',14);
        colormap(jet);
        
        if SAVE
            saveas(gcf,[filepath '/' name '_uc' num2str(updist*1e6)  '.png'])
            if exist('newLED','var')
                if name == 'B111dataToPlot'
                    B111ferro = B_data;
                    eval(['save ' filepath filesep name '_uc' num2str(updist*1e6) '.mat B111ferro Bt h step newLED corners -mat']);
                else
                    Bz = B_data;
                    eval(['save ' filepath filesep name '_uc' num2str(updist*1e6) '.mat Bz Bt h step newLED corners -mat']);
                end
            else
                if name == 'B111dataToPlot'
                    B111ferro = B_data;
                    eval(['save ' filepath filesep name '_uc' num2str(updist*1e6) '.mat B111ferro h step -mat']);
                else
                    Bz = B_data;
                    eval(['save ' filepath filesep name '_uc' num2str(updist*1e6) '.mat Bz h step -mat']);
                end
            end
            disp(sprintf(['\nSaving' name '_uc' num2str(updist*1e6) '.mat...\n']))
        end
    end
end
