function  [resid, Bzmodel, M] = SourceFitMultiP8(P,x,y,BzExp,display,method,quad,morder);
%{

parameters:
    P: double
        initial guess, with [x0, y0, h]
    x: []
        X-coordinates
    y: []
        Y-coordinates
    BzExp:
        Bz data
    display: bool
        if 1 plots Bz and residual
    method: bool
        if true uses norm 
    quad: quadrupole?
    morder: 1, 2 or 3
        if 1: dipole fit
        if 2: quadrupole fit
        if 3: octopole fit
%}
COIL=0;             %1 = circular coil, 0 = square coil
THREED=0;
MINNORM=0;
TERMS=[3 8 15];
const = 1.0*10^-7; %field in T

Rcoil=40e-6; %75e-6; %120e-6;


x0=P(1);
y0=P(2);
h=P(3);

% if quad
%      Bz=zeros(size(BzExp));
%      if COIL
%          w1=0.785398163397449;                   %9-point quadrature on a disk with radius Rcoil
%          w2=0.294524311274043;
%          D1=0.754344479484572;
%          D2=0.312459714103782;
%          wi=[w1 w2 w2 w2 w2 w2 w2 w2 w2]/pi;
%          DLx=[0 -D1 -D1  D1  D1 -D2 -D2  D2  D2]*Rcoil*0.5;
%          DLy=[0 -D2  D2 -D2  D2 -D1  D1 -D1  D1]*Rcoil;
%          
%          for k1=1:9
%              delx = (x-DLx(k1)) - x0;
%              dely = (y-DLy(k1)) - y0;
%              l5i = (delx.^2+dely.^2+h^2).^-2.5;
%              l3i = (delx.^2+dely.^2+h^2).^-1.5;
%              Bz=Bz+wi(k1)*const*3*h*m*((delx*sind(theta)*cosd(phi)+ dely*sind(theta)*sind(phi) + h*cosd(theta)).*l5i-cosd(theta).*l3i/3/h); %field in nT
%          end
%      else
%          wix=[5/9 8/9 5/9]/2;                     %9-point quadrature on a square of size 2*Rcoil
%          wiy=[5/9 8/9 5/9]/2;
%          DLx=[-sqrt(3/5) 0 sqrt(3/5)]*Rcoil*0.5;
%          DLy=[-sqrt(3/5) 0 sqrt(3/5)]*Rcoil;
%          h0=h;
%          for k1=1:3         %top layer of the cube
%              for k2=1:3
%                  delx = (x-DLx(k1)) - x0;
%                  dely = (y-DLy(k2)) - y0;
%                  l5i = (delx.^2+dely.^2+h^2).^-2.5;
%                  l3i = (delx.^2+dely.^2+h^2).^-1.5;
%                  %Bz=Bz+wix(k1)*wiy(k2);
%                  Bz=Bz+wix(k1)*wiy(k2)*const*3*h*m*((delx*sind(theta)*cosd(phi)+ dely*sind(theta)*sind(phi) + h*cosd(theta)).*l5i-cosd(theta).*l3i/3/h); %field in nT
%              end
%          end
%          if THREED
%              h=h+Rcoil/2;         
%              for k1=1:3         %middle layer of the cube
%                  for k2=1:3
%                      delx = (x-DLx(k1)) - x0;
%                      dely = (y-DLy(k2)) - y0;
%                      l5i = (delx.^2+dely.^2+h^2).^-2.5;
%                      l3i = (delx.^2+dely.^2+h^2).^-1.5;
%                      %Bz=Bz+wix(k1)*wiy(k2);
%                      Bz=Bz+wix(k1)*wiy(k2)*const*3*h*m*((delx*sind(theta)*cosd(phi)+ dely*sind(theta)*sind(phi) + h*cosd(theta)).*l5i-cosd(theta).*l3i/3/h); %field in nT
%                  end
%              end
%              h=h+Rcoil/2;
%              for k1=1:3         %bottom layer of the cube
%                  for k2=1:3
%                      delx = (x-DLx(k1)) - x0;
%                      dely = (y-DLy(k2)) - y0;
%                      l5i = (delx.^2+dely.^2+h^2).^-2.5;
%                      l3i = (delx.^2+dely.^2+h^2).^-1.5;
%                      %Bz=Bz+wix(k1)*wiy(k2);
%                      Bz=Bz+wix(k1)*wiy(k2)*const*3*h*m*((delx*sind(theta)*cosd(phi)+ dely*sind(theta)*sind(phi) + h*cosd(theta)).*l5i-cosd(theta).*l3i/3/h); %field in nT
%                  end
%              end
%              Bz=Bz/3; %normalize by the number of layers
%              h=h0;
%          end
%      end
% else

if quad && COIL
    w1=0.785398163397449;                   %9-point quadrature on a disk with radius Rcoil
    w2=0.294524311274043;
    D1=0.754344479484572;
    D2=0.312459714103782;
    wi=[w1 w2 w2 w2 w2 w2 w2 w2 w2]/pi;
    DLx=[0 -D1 -D1  D1  D1 -D2 -D2  D2  D2]*Rcoil;
    DLy=[0 -D2  D2 -D2  D2 -D1  D1 -D1  D1]*Rcoil;
    A=zeros(numel(BzExp),15);
    for k1=1:9
        delx = (x-DLx(k1)) - x0;
        dely = (y-DLy(k1)) - y0;
        l5i = (delx.^2+dely.^2+h^2).^-2.5;
        %dipole terms
        A(:,1)=A(:,1)+wi(k1)*3*const*h*delx(:).*l5i(:);
        A(:,2)=A(:,2)+wi(k1)*3*const*h*dely(:).*l5i(:);
        A(:,3)=A(:,3)+wi(k1)*const*(2*h*h*ones(numel(BzExp),1)-delx(:).^2-dely(:).^2).*l5i(:);
        %quadrupole terms
        l7i = (delx.^2+dely.^2+h^2).^-3.5;
        A(:,4)=A(:,4)-wi(k1)*const*( (9*h).*((1/2)*l5i(:)) - (15*h^3).*((1/2)*l7i(:)) );   %Phi_e,20
        A(:,5)=A(:,5)-wi(k1)*const*( (3*delx(:)).*l5i(:) - (15*delx(:).*h^2).*l7i(:) );   %Phi_e,21
        A(:,6)=A(:,6)-wi(k1)*const*( (3*dely(:)).*l5i(:) - (15*dely(:).*h^2).*l7i(:) );   %Phi_o,21
        A(:,7)=A(:,7)-wi(k1)*const*( -(5*h*(3*delx(:).^2 - 3*dely(:).^2)).*l7i(:) );   %Phi_e,22
        A(:,8)=A(:,8)-wi(k1)*const*( -(30*delx(:).*dely(:)*h).*l7i(:) );   %Phi_o,22
        
        %octupole terms
        l9i = (delx.^2+dely.^2+h^2).^-4.5;
        A(:,9)=A(:,9)+wi(k1)*const*((3*delx(:).^4 + 6*delx(:).^2.*dely(:).^2 - 24*delx(:).^2*h^2 + 3*dely(:).^4 - 24*dely(:).^2*h^2 + 8*h^4).*l9i(:))/2;
        A(:,10)=A(:,10)-wi(k1)*const*((15*delx(:)*h.*(3*delx(:).^2 + 3*dely(:).^2 - 4*h^2)).*l9i(:))/2;
        A(:,11)=A(:,11)-wi(k1)*const*((15*dely(:)*h.*(3*delx(:).^2 + 3*dely(:).^2 - 4*h^2)).*l9i(:))/2;
        A(:,12)=A(:,12)-wi(k1)*const*((15*(delx(:) - dely(:)).*((delx(:) + dely(:)).*(delx(:).^2 + dely(:).^2 - 6*h^2))).*l9i(:));
        A(:,13)=A(:,13)-wi(k1)*const*((30*delx(:).*dely(:).*(delx(:).^2 + dely(:).^2 - 6*h^2)).*l9i(:));
        A(:,14)=A(:,14)+wi(k1)*const*((105*delx(:)*h.*(delx(:).^2 - 3*dely(:).^2)).*l9i(:));
        A(:,15)=A(:,15)+wi(k1)*const*((105*dely(:)*h.*(3*delx(:).^2 - dely(:).^2)).*l9i(:));
    end
else
    A=zeros(numel(BzExp),TERMS(morder));
    delx = x-x0;
    dely = y-y0;
    l5i = (delx.^2+dely.^2+h^2).^-2.5;
    %dipole terms
    A(:,1)=3*const*h*delx(:).*l5i(:);
    A(:,2)=3*const*h*dely(:).*l5i(:);
    A(:,3)=const*(2*h*h*ones(numel(BzExp),1)-delx(:).^2-dely(:).^2).*l5i(:);
    
    if morder>1
        %quadrupole terms
        l7i = (delx.^2+dely.^2+h^2).^-3.5;
        A(:,4)=-const*( (9*h).*((1/2)*l5i(:)) - (15*h^3).*((1/2)*l7i(:)) );   %Phi_e,20
        A(:,5)=-const*( (3*delx(:)).*l5i(:) - (15*delx(:).*h^2).*l7i(:) );   %Phi_e,21
        A(:,6)=-const*( (3*dely(:)).*l5i(:) - (15*dely(:).*h^2).*l7i(:) );   %Phi_o,21
        A(:,7)=-const*( -(5*h*(3*delx(:).^2 - 3*dely(:).^2)).*l7i(:) );   %Phi_e,22
        A(:,8)=-const*( -(30*delx(:).*dely(:)*h).*l7i(:) );   %Phi_o,22
    end
    
    if morder>2
        %octupole terms
        l9i = (delx.^2+dely.^2+h^2).^-4.5;
        A(:,9)=const*((3*delx(:).^4 + 6*delx(:).^2.*dely(:).^2 - 24*delx(:).^2*h^2 + 3*dely(:).^4 - 24*dely(:).^2*h^2 + 8*h^4).*l9i(:))/2;
        A(:,10)=-const*((15*delx(:)*h.*(3*delx(:).^2 + 3*dely(:).^2 - 4*h^2)).*l9i(:))/2;
        A(:,11)=-const*((15*dely(:)*h.*(3*delx(:).^2 + 3*dely(:).^2 - 4*h^2)).*l9i(:))/2;
        A(:,12)=-const*((15*(delx(:) - dely(:)).*((delx(:) + dely(:)).*(delx(:).^2 + dely(:).^2 - 6*h^2))).*l9i(:));
        A(:,13)=-const*((30*delx(:).*dely(:).*(delx(:).^2 + dely(:).^2 - 6*h^2)).*l9i(:));
        A(:,14)=const*((105*delx(:)*h.*(delx(:).^2 - 3*dely(:).^2)).*l9i(:));
        A(:,15)=const*((105*dely(:)*h.*(3*delx(:).^2 - dely(:).^2)).*l9i(:));
    end
end

    %    warning off
    %keyboard
    Scl=diag(1./sqrt(sum(A.^2)));
    Ascl=A*Scl;
    Mu=Ascl\BzExp(:);
    
    M=Scl*Mu;
%    warning on 
 
    Bz=A*M; %field in nT
    Bz=reshape(Bz,size(BzExp));
    %keyboard
%end

residual=Bz-BzExp;

if display
    sourceFitFigure = figure();
    subplot(2,2,2);
    imagesc(Bz);
    axis xy, axis equal, axis tight
    caxis([-1 1]*max(abs(caxis)));
    colorbar
    subplot(2,2,3);
    imagesc(residual);
    axis xy, axis equal, axis tight
    caxis([-1 1]*max(abs(caxis)));
    colorbar
    drawnow;
    %pause(1)
end

if method
    if MINNORM
        resid=norm(residual,2)+gamma*norm(Mu,2);
    else
        resid=norm(residual,2);
    end
else
    if MINNORM
        resid=residual*(norm(Mu,2)^(1/2));
    else
        resid=residual;
    end
end
%resid=norm(residual,2);
Bzmodel=Bz;
end
