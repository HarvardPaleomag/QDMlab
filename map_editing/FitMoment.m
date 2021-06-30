function residex=FitMoment(MORDER,INFILE, XY1, XY2, NRUNS, QUAD)
%MORDER is highest order in the multipole expansion: 2 is quadrupole; 3 is
%octupole
%INFILE is the name to the .mat data file with Bz map
%XY1 and XY2 are in form [x1 y1] and [x2 y2]

rng('shuffle')

set(0,'DefaultFigureColormap',jet)

OutFileName='SubtractedSources.txt';

%--------------------------------------------------------------------------

%Initial parameters guess (before randomization)

incl=-0;
dec=0;
%incl=40;
%dec=300;

m0=1e-12;

METHOD=0;       %0 = least squares, 1=Nelder-Mead
SNR=0;          %signal-to-noise ratio in dB
TERMS=[3,8,15];
MINTOL=1;
NSTAT=50;
DISPLAY=0;

theta0=90-incl;
phi0=dec+90;    %SM y axis is the one pointing northy


%----------------------------------------------------------------------


counter=1;

while counter
    
    step=0;
    load(INFILE);
    
    %crop maps to the specified region
    Bz=Bz(XY1(1,2):XY2(1,2),XY1(1,1):XY2(1,1));
    Bt=Bt(XY1(1,2):XY2(1,2),XY1(1,1):XY2(1,1));
    
    XY1
    XY2
    
    %note that "scanc" already holds a cropped map
    scanc=Bz;
    
    %make array of x and y with physical units
    xc=((1:size(scanc,2)))*step;
    yc=((1:size(scanc,1)))*step;
    [Xc,Yc]=meshgrid(xc,yc);
    
    P00(1)=round(size(scanc,2)/2)*step;
    P00(2)=round(size(scanc,1)/2)*step;
    P00(3)=h;
    
    P=zeros(length(P00)+TERMS(MORDER),NRUNS);
    fval=zeros(1,NRUNS);
    fval2=zeros(1,NRUNS);
    for k=1:NRUNS
        if NRUNS==1
            P0=P00;
        else
            P0=P00+0.2*(rand(size(P00))-0.5).*P00;
        end
        
        options=optimset('TolX',10^(floor(log10(m0))-5),'TolFun',10^(floor(log10(max(abs(Bz(:)))))-8),'MaxFunEvals',6000,'MaxIter',2000,'Display','none');
        
        if METHOD
            [P(1:3,k),fval2(k),exitflag,output] = fminsearch(@(Pp) SourceFitMultiP8(Pp,Xc,Yc,scanc,DISPLAY,METHOD,QUAD,MORDER),P0,options);
        else
            [P(1:3,k),fval2(k),resd,exitflag,output] = lsqnonlin(@(Pp) SourceFitMultiP8(Pp,Xc,Yc,scanc,DISPLAY,METHOD,QUAD,MORDER),P0,[],[],options);
        end
        [resid,BzModel,M]=SourceFitMultiP8(P(1:3,k),Xc,Yc,scanc,DISPLAY,METHOD,QUAD,MORDER);
        
        Mx=M(1);
        My=M(2);
        Mz=M(3);
        m=sqrt(Mx^2+My^2+Mz^2);
        theta=acosd(Mz/m);
        phi=atan2d(My,Mx);
        P(4,k)=m;
        %convert angles to inclination and declination
        P(5,k)=90-theta;
        P(6,k)=phi-90;
        
        %enforce range for parameters
        change=1;
        while change
            change=0;
            if P(5,k)>90
                P(5,k)=180-P(5,k); % i' = 180-i
                P(6,k)=P(6,k)+180; % d' = d+180
                change=1;
            end
            if P(5,k)<-90
                P(5,k)=-180-P(5,k); % i' = -180-i
                P(6,k)=P(6,k)+180; % d' = d+180
                change=1;
            end
            if P(6,k)<0
                P(6,k)=P(6,k)+360; % d' = d+360
                change=1;
            end
            if P(6,k)>=360
                P(6,k)=P(6,k)-360; % d' = d-360
                change=1;
            end
        end
        if MORDER>1
            %Quadrupole moment
            P(7,k)=M(4);
            P(8,k)=M(5);
            P(9,k)=M(6);
            P(10,k)=M(7);
            P(11,k)=M(8);
        end
        
        if MORDER>2
            %Octupole moment
            P(12,k)=M(9);
            P(13,k)=M(10);
            P(14,k)=M(11);
            P(15,k)=M(12);
            P(16,k)=M(13);
            P(17,k)=M(14);
            P(18,k)=M(15);
        end
        
        fprintf('..%0d..(%0d)  ',k,output.iterations);
        if ~mod(k,10)
            fprintf('\n');
        end
        if size(P,1)>6
            Paux=[P(1:4,k)' 90-P(5,k) P(6,k)+90 P(7:TERMS(MORDER)+3,k)'];
        else
            Paux=[P(1:4,k)' 90-P(5,k) P(6,k)+90];
        end
        [resid,Bzmodel]=SourceFitMultiP8(Paux,Xc,Yc,scanc,0,METHOD,QUAD,MORDER);
        fval(k)=sqrt(sum(sum((Bzmodel-scanc).^2))/numel(scanc));
    end
    i0=find(fval==min(fval));
    fsort=sort(fval);
    
    %i=find(fval<=min(fval)+MINTOL*(max(fval)-min(fval)));
    i=find(fval<=fsort(MINTOL));
    disp(sprintf('\n--- Averaging %d points ---',numel(i)))
    if numel(i)>0.1*numel(fval)
        disp('Too many points are being averaged. Consider adjusting MINTOL parameter.')
    end
    
    
    Popt=zeros(size(Paux));
    Popt(4)=sum(P(4,i).*fval(i))/sum(fval(i));
    mopt=Popt(4)
    
    Popt(5)=sum(P(5,i).*fval(i))/sum(fval(i));
    iopt=Popt(5);
    
    Popt(6)=sum(P(6,i).*fval(i))/sum(fval(i));
    dopt=Popt(6);
    
    Popt(1)=sum(P(1,i).*fval(i))/sum(fval(i));
    Popt(2)=sum(P(2,i).*fval(i))/sum(fval(i));
    Popt(3)=sum(P(3,i).*fval(i))/sum(fval(i));
    hopt=Popt(3);
    
    for kk=7:TERMS(MORDER)+3
        Popt(kk)=sum(P(kk,i).*fval(i))/sum(fval(i));
    end
    
    xopt=Popt(1);
    
    yopt=Popt(2);
    
    Popt2=[Popt(1:4) 90-Popt(5) Popt(6)+90 Popt(7:TERMS(MORDER)+3)];
    [resid,Bzmodel]=SourceFitMultiP8(Popt2,Xc,Yc,scanc,0,METHOD,QUAD,MORDER);
    residex=Bzmodel-scanc;
    
    figure
    subplot(2,2,1);
    imagesc(xc,yc,scanc);
    axis xy, axis equal, axis tight;
    caxis([-1 1]*max(abs(caxis)));
    colorbar
    title('Original Scan');
    subplot(2,2,2);
    imagesc(xc,yc,Bzmodel);
    axis xy, axis equal, axis tight;
    caxis([-1 1]*max(abs(caxis)));
    colorbar
    title('Model Scan');
    subplot(2,2,4);
    imagesc(xc,yc,residex);
    axis xy, axis equal, axis tight;
    caxis([-1 1]*max(abs(caxis)));
    colorbar
    title('Residuals');
    
    %resids=sqrt(sum(sum(residex.^2))/numel(residex))
    resids=sqrt(sum(sum(residex.^2))/sum(sum(scanc.^2)));
    [path,name,ext]=fileparts(INFILE);
    nameext=[name ext];
    
    fid=fopen([path '/' OutFileName],'r');
    header=(fid==-1);
    if ~header
        fclose(fid);
    end
    fid=fopen([path '/' OutFileName],'a+t');
    if header
        fprintf(fid,'File Name\tMoment\tInclination\tDeclination\tHeight\tResiduals\r\n');
    end
    %Note a 180 rotation about y axis is imposed here
    if hopt < 0
        fprintf(fid,'%s\t%1.5d\t%1.5d\t%1.5d\t%1.5d\t%1.5d\r\n',nameext,mopt,-iopt,mod(180-dopt,360),-hopt,resids);
    else
        fprintf(fid,'%s\t%1.5d\t%1.5d\t%1.5d\t%1.5d\t%1.5d\r\n',nameext,mopt,-iopt,mod(360-dopt,360),hopt,resids);
    end
    fclose(fid);
    
    counter=counter-1;
    
end


