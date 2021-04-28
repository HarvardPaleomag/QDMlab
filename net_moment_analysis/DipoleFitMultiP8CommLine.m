function results = DipoleFitMultiP8CommLine(MORDER,INFILE, XY, CROPFACT, DOWNSMPL, NRUNS, QUAD, outputtrue,plottrue, varargin)
%MORDER is highest order in the multipole expansion: 2 is quadrupole; 8 is
%octupole
% 
% 
% Parameters
% ----------
%     positional
%     ----------
%     MORDER:
%     INFILE:
%     XY:
%     CROPFACT:
%     DOWNSMPL:
%     NRUNS:
%     QUAD:
%     outputtrue:
%     plottrue:
%     
%     optional parameters
%     -------------------
%     dx: int
%         default: 0
%         the width of the box that is cropped around the data
%     dy: int
%         default: 0
%         the height of the box that is cropped around the data
%     data: struct
%         default: 0
%         If you want to load data directly instead of reloading it from the
%         file you can pass it after the 'data' keyword.
%     
%     following parameters for costraining fit
%     ========================================
%     m0:
%         default: 1e-12
%     hguess:
%         default: 5e-5
%     minheight:
%         default: 0
%     maxheight:
%         default: 100e-6
%     boxwidth:
%         default: 100e-6
% 
% 
%     STATISTSICS:
%         default: 0
%     NSTAT:
%         default: 50
%     DISPLAY:
%         default: 0
% 
% Returns: struct
%     returns a structure with the following entries:
%         dfile: file path
%         m: fitted moment
%         inc: fitted inclination
%         dec: fitted declination
%         h: fitted height (depth) of dipole
%         res: residuals


% add optional arguments to the function by using inputParser
inParams = inputParser;% define parameters
addRequired(inParams,'MORDER');
addRequired(inParams,'INFILE');
addRequired(inParams,'XY');
addRequired(inParams,'CROPFACT');
addRequired(inParams, 'DOWNSMPL');
addRequired(inParams, 'NRUNS');
addRequired(inParams, 'QUAD');
addRequired(inParams, 'outputtrue');
addRequired(inParams, 'plottrue');
addParameter(inParams, 'dx', 0)
addParameter(inParams, 'dy', 0)
addParameter(inParams, 'data', 0)
addParameter(inParams, 'm0', 1e-12)
addParameter(inParams, 'hguess', 5e-5)
addParameter(inParams, 'minheight', 0)
addParameter(inParams, 'maxheight', 100e-6)
addParameter(inParams, 'boxwidth', 100e-6)
addParameter(inParams, 'STATISTICS', 0)
addParameter(inParams, 'NSTAT', 50)
addParameter(inParams, 'DISPLAY', 0)
addParameter(inParams, 'IMAGEFOLDER',0)
addParameter(inParams, 'SOURCENAME',0)

parse(inParams,MORDER,INFILE, XY, CROPFACT, DOWNSMPL, NRUNS, QUAD, outputtrue,plottrue,varargin{:});

% create variables from the optional parameters of the parser
dx = inParams.Results.dx;
dy = inParams.Results.dy;
data = inParams.Results.data;
m0=inParams.Results.m0;
hguess=inParams.Results.hguess;
minheight=inParams.Results.minheight;
maxheight=inParams.Results.maxheight;
boxwidth=inParams.Results.boxwidth;
STATISTICS=inParams.Results.STATISTICS;
NSTAT=inParams.Results.NSTAT;
DISPLAY=inParams.Results.DISPLAY;
IMAGEFOLDER=inParams.Results.IMAGEFOLDER;
SOURCENAME=inParams.Results.SOURCENAME;

rng('shuffle')

set(0,'DefaultFigureColormap',jet)

OutFileName='DipoleInversions.txt';

%--------------------------------------------------------------------------

%Initial parameters guess (before randomization)

incl=-0;
dec=0;
%incl=40;
%dec=300;



METHOD=0;       %0 = least squares, 1=Nelder-Mead
NOISE=0;        %0 = no noise
SNR=0;          %signal-to-noise ratio in dB
AUTO=0;         %automatically find dipole position from Bt map
TERMS=[3,8,15];
MINTOL=1;

theta0=90-incl;
phi0=dec+90;    %SM y axis is the one pointing northy


%----------------------------------------------------------------------

if STATISTICS
    counter=NSTAT;
else
    counter=1;
end

%%
while counter

    if STATISTICS
        disp(sprintf('* Iteration %d...',NSTAT-counter+1));
    end
    step=0;

    % check if the optional parameter 'data' was passed.
    % if yes, newLED, Bz, h, and step variables need to be assigned
    if isa(data, 'struct')
        newLED = data.newLED;
        Bz = data.Bz;
        
        h = data.h;
        step = data.step; % size of a pixel?
        
        % check if h/step is too large
        if h > 0.1
            h = h*1e-6;
        end
        
        if step > 1e-2
            step = step*1e-6;
        end

    else
        disp('<> LOADING.. file')
        load(INFILE);
    end
    
    [filepath,name,ext]=fileparts(INFILE);
   
    if DOWNSMPL
        Bz=downsample(downsample(Bz,DOWNSMPL)',DOWNSMPL)';
        step=step*DOWNSMPL;
    end
    
    scan=Bz;    % Bz is assumed in T
    h=h;%*1e-6;
    disp(h)
    x=((1:size(scan,2))-1)*step;
    y=((1:size(scan,1))-1)*step;
    [X,Y]=meshgrid(x,y);
    
    if ~STATISTICS
        figure(1);
    end
    
    if NOISE
        %noise=randn(size(scan))*0.5*std(std(scan));
        noise=sqrt(10^(-SNR/10)*var(scan(:)))*randn(size(scan));
        scan=scan+noise;
    end
    
    if DISPLAY
        NRUNS=1;
    end

    
    % i, j are the approximate center of the dipole in map coordinates
    j = round(XY(1)/DOWNSMPL);
    i = round(XY(2)/DOWNSMPL);
    
    % cropi, cropj are the cropping indices for x(i) and y(j) of the
    % original map
    
    % check if the optional parameters dx and dy are not 0
    if all([dx,dy])
        % define cropj, cropi from dx, dy instead of CROPFACT
        cropj = round([j j+dx]);
        cropi = round([i i+dy]);
        i = i+dy/2; j = j+dx/2;
    else
        cropj = round(1+j+[-CROPFACT CROPFACT]);
        cropi = round(1+i+[-CROPFACT CROPFACT]);
    end
    
    % position of dipole
    x0 = j*step;
    y0 = i*step;
    
    %% adjust if the crop area falls outside the image
    scansize = size(scan);
    for p=1:2
        if cropj(p) < 1
            cropj(p) = 1;
        end
        if cropi(p) < 1
            cropi(p) = 1;
        end
        if cropj(p) > scansize(2)
            cropj(p) = scansize(2);
        end
        if cropi(p) > scansize(1)
            cropi(p) = scansize(1);
        end
    end
    
    scanCropped = scan(cropi(1):cropi(2),cropj(1):cropj(2));
    Xc=X(cropi(1):cropi(2),cropj(1):cropj(2));
    Yc=Y(cropi(1):cropi(2),cropj(1):cropj(2));
    xc=x(cropj(1):cropj(2));
    yc=y(cropi(1):cropi(2));
    
    figure
    imagesc(xc,yc,abs(Bz(cropi(1):cropi(2),cropj(1):cropj(2))));
    axis xy, axis equal, axis tight
    caxis([0 1]*max(abs(caxis)));
    colormap(hot)
    colorbar
    title(sprintf('Cropped map'));
    drawnow
    
    P00(1)=x0;
    P00(2)=y0;
    P00(3)=hguess;
    
    drawnow
    P=zeros(length(P00)+TERMS(MORDER),NRUNS);
    fval=zeros(1,NRUNS);
    fval2=zeros(1,NRUNS);
    
    for k=1:NRUNS
        if NRUNS==1
            P0=P00;%+0.3*(rand(size(P00))-0.5).*P00;
        else
            %P0=P00+0.1*(rand(size(P00))-0.5).*[P00(1) 20 40 P00(4:6)];
            P0=P00+0.1*(rand(size(P00))-0.5).*P00;
        end
        options=optimset('TolX',10^(floor(log10(m0))-5),'TolFun', 10^(floor(log10(max(abs(Bz(:)))))-8),'MaxFunEvals',6000,'MaxIter',2000,'Display','none');
        
        if METHOD
            [P(1:3,k),fval2(k),exitflag,output] = fminsearch(@(Pp) SourceFitMultiP8(Pp,Xc,Yc,scanCropped,DISPLAY,METHOD,QUAD,MORDER),P0,options);

%             [P(1:3,k),fval2(k),exitflag,output] = fmincon(@(Pp) SourceFitMultiP8(Pp,Xc,Yc,scanCropped,DISPLAY,METHOD,QUAD,MORDER),P0,[],[],[],[],[x0-boxwidth,y0-boxwidth,minheight],[x0+boxwidth,y0+boxwidth,maxheight],[],options);
        else
            [P(1:3,k),fval2(k),resd,exitflag,output] = lsqnonlin(@(Pp) SourceFitMultiP8(Pp,Xc,Yc,scanCropped,DISPLAY,METHOD,QUAD,MORDER),P0,[],[],options);
        end
        [resid,BzModel,M]=SourceFitMultiP8(P(1:3,k),Xc,Yc,scanCropped,DISPLAY,METHOD,QUAD,MORDER);
        
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
        if DISPLAY
            disp(sprintf('Moment: %.4e',P(4)))
            disp(sprintf('Inclination: %.1f',P(5)))
            disp(sprintf('Declination: %.1f',P(6)))
            disp(sprintf('Height: %.4e',P(3)))
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
        [resid,Bzmodel]=SourceFitMultiP8(Paux,Xc,Yc,scanCropped,0,METHOD,QUAD,MORDER);
        % calculate RMS of model-scan
        fval(k)=sqrt(sum(sum((Bzmodel-scanCropped).^2))/numel(scanCropped));
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
    mopt=Popt(4)%*1000
    disp(sprintf('(min = %1.3d)',P(4,i0)));
    
    Popt(5)=sum(P(5,i).*fval(i))/sum(fval(i));
    iopt=Popt(5);
    -iopt
    disp(sprintf('(min = %1.3f)',P(5,i0)));
    
    Popt(6)=sum(P(6,i).*fval(i))/sum(fval(i));
    dopt=Popt(6);                 %x and y axis are reversed in SM - N is Y instead of x
    mod(360-dopt,360)
    disp(sprintf('(min = %1.3f)',P(6,i0)));
    
    Popt(1)=sum(P(1,i).*fval(i))/sum(fval(i));
    Popt(2)=sum(P(2,i).*fval(i))/sum(fval(i));
    Popt(3)=sum(P(3,i).*fval(i))/sum(fval(i));
    hopt=Popt(3)
    disp(sprintf('(min = %1.3d)',P(3,i0)));
    
    for kk=7:TERMS(MORDER)+3
        Popt(kk)=sum(P(kk,i).*fval(i))/sum(fval(i));
    end    
    
    xopt=Popt(1)
    disp(sprintf('(min = %1.3d)',P(1,i0)));
    
    yopt=Popt(2)
    disp(sprintf('(min = %1.3d)',P(2,i0)));
    
    if plottrue
        if ~STATISTICS
            figure
            plot(P(1,:),fval,'.')
            title('Moment');
            hold on
            plot(P(1,i),fval(i),'r.')
            plot(P(1,i0),fval(i0),'c.')
            plot([mopt mopt],ylim,'m--');
            hold off
            
            figure
            plot(P(2,:),fval,'.')
            title('Inclination');
            hold on
            plot(P(2,i),fval(i),'r.')
            plot(P(2,i0),fval(i0),'c.')
            plot([iopt iopt],ylim,'m--');
            hold off
            
            figure
            plot(P(3,:),fval,'.')
            title('Declination');
            hold on
            plot(P(3,i),fval(i),'r.')
            plot(P(3,i0),fval(i0),'c.')
            plot([dopt dopt],ylim,'m--');
            hold off
            
            figure
            plot(P(6,:),fval,'.')
            title('Height');
            hold on
            plot(P(6,i),fval(i),'r.')
            plot(P(6,i0),fval(i0),'c.')
            plot([hopt hopt],ylim,'m--');
            hold off
            
            figure
            plot(P(4,:),fval,'.')
            title('X displacement');
            hold on
            plot(P(4,i),fval(i),'r.')
            plot(P(4,i0),fval(i0),'c.')
            plot([xopt xopt],ylim,'m--');
            hold off
            
            figure
            plot(P(5,:),fval,'.')
            title('Y displacement');
            hold on
            plot(P(5,i),fval(i),'r.')
            plot(P(5,i0),fval(i0),'c.')
            plot([yopt yopt],ylim,'m--');
            hold off
        end
    end
    
    Popt2=[Popt(1:4) 90-Popt(5) Popt(6)+90 Popt(7:TERMS(MORDER)+3)];
    [resid,Bzmodel]=SourceFitMultiP8(Popt2,Xc,Yc,scanCropped,0,METHOD,QUAD,MORDER);
    residex=Bzmodel-scanCropped;
    if ~STATISTICS
        figure
        subplot(2,2,1);
        imagesc(xc,yc,scanCropped);
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
    end
    if IMAGEFOLDER
        if SOURCENAME
            saveas(gcf,[IMAGEFOLDER '/Fit_' name '_M' num2str(MORDER) '_' SOURCENAME '.png'])
        else
            disp('Saving to IMAGEFOLDER directory requires a SOURCENAME parameter.');
        end
    else
        saveas(gcf,[filepath '/Fit_' name '_M' num2str(MORDER) '_x' num2str(round(XY(1))) 'y' num2str(round(XY(2))) '.png'])
    end
    
    %resids=sqrt(sum(sum(residex.^2))/numel(residex))
    resids=sqrt(sum(sum(residex.^2))/sum(sum(scanCropped.^2)));
    [path,name,ext]=fileparts(INFILE);
    nameext=[name ext];
    
    if STATISTICS
        outputtrue=1;
    end
   
    if outputtrue
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
    end
    
    counter=counter-1;
    
end

% create the outputs of the funtion
if hopt < 0
    results = struct('dfile',INFILE, 'm', mopt, 'inc', -iopt, 'dec', mod(180-dopt,360), 'h', -hopt, 'res', resids);
else
    results = struct('dfile',INFILE, 'm', mopt, 'inc', -iopt, 'dec', mod(360-dopt,360), 'h', -hopt, 'res', resids);
end

