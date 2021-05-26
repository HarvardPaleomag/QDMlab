function [Bout] = UpCont(Bin, dh, fs)
%[Bout] = UpCont(Bin, dh, fs)
% Parameters
% ----------
%       Bin: double
%           component map
%       dh: double
%           distance to upward continue the field map
%       fs: double
%           Sampling frequency in 1/m

EXPAND = 1;
FOVERSAMPL = 2; % Algorithm refines sampling in the frequency domain by
% this factor (i.e., zero padding). Increase value (should be
if EXPAND
    Bino = Bin; % an even number) if the number of points in Bz map is small.
    % bkg=Bin(:,1)*ones(1,size(Bin,2));
    % Bin=Bin-bkg;
    % Bin=Bin';
    Bin = [zeros(size(Bin)), zeros(size(Bin)), zeros(size(Bin)); ...
        zeros(size(Bin)), Bin, zeros(size(Bin)); ...
        zeros(size(Bin)), zeros(size(Bin)), zeros(size(Bin))];
end

[SIZEx, SIZEy] = size(Bin');
Nx = SIZEx * FOVERSAMPL;
Ny = SIZEy * FOVERSAMPL;
%fs= fs * FOVERSAMPL;                        % Adjust sampling frequency accordingly


f2 = [0:Ny / 2, -(Ny / 2 - 1):-1] * fs / Ny; % These freq. coordinates match MATLAB's fft algorithm
f1 = [0:Nx / 2, -(Nx / 2 - 1):-1] * fs / Nx;

[F2, F1] = meshgrid(f2, f1); % Generate spatial frequency grid. A tiny 'epsilon' is added to
% each freq. variable so as to avoid division by zero warning.
kx = 2 * pi * F1;
ky = 2 * pi * F2;
k = sqrt(kx.^2+ky.^2);

et = exp(-dh*k); % Calculate the filter frequency response associated with the x component


e = fft2(Bin', Nx, Ny); % Compute FFT of the field map

Bout = ifft2(e.*et, Nx, Ny, 'symmetric')'; % Calculate x component


Bout = Bout(1:Ny/FOVERSAMPL, 1:Nx/FOVERSAMPL); % Crop matrices to get rid of zero padding
%Bout=Bout-mean2(Bout);
if EXPAND
    Bout = Bout(size(Bino, 1)+1:2*size(Bino, 1), size(Bino, 2)+1:2*size(Bino, 2));
end
