function [Bz] = QDMBzFromBu(Bu, fs, u)
%[Bz] = QDMBzFromBu(Bu, fs, u)
%   Retrieves the z component of the magnetic field from a QDM map of the u component. It performs
%   this operation in the frequency domain.
%
%   Code by Eduardo A. Lima - (c) 2017

%   ----------------------------------------------------------------------------------------------
%   Bu     -> u component of the magnetic field measured in a regular planar grid
%   fs     -> Sampling frequency in 1/m
%   u      -> unit vector representing the field component measured
%   ----------------------------------------------------------------------------------------------


SHOWGRAPHS = 0; % Set this value to 1 to show frequency response plots, or to 0 otherwise.

[SIZEx, SIZEy] = size(Bu);
N1 = SIZEx;
N2 = SIZEy;

f1 = [0:N1 / 2 - 1, -(N1 / 2):-1] * fs / N1; %these freq. coordinates match the fft algorithm
f2 = [0:N2 / 2 - 1, -(N2 / 2):-1] * fs / N2;

[F2, F1] = meshgrid(f2+1e-30, f1+1e-30);

ky = 2 * pi * F1;
kx = 2 * pi * F2;

k = sqrt(kx.^2+ky.^2);


etz = k ./ (u(3) * k - u(2) * 1i * ky - u(1) * 1i * kx); % calculate the filter frequency response associated with the x component


e = fft2(Bu, N1, N2);

Bz = ifft2(e.*etz, N1, N2, 'symmetric');
%Bz=real(ifft2(e.*etz));
