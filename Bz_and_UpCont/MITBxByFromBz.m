function [Bx, By] = MITBxByFromBz(Bz, fs)
%[Bx, By] = MITBxByFromBz(Bz, fs)
%   Retrieves the x and y components of the magnetic field from a map of the z component. It performs
%   this operation in the frequency domain. See, for instance, Egli and Heller (eq. 12), and Roth and
%   Wikswo.
%
%   Code by Eduardo A. Lima - (c) 2007

%   ----------------------------------------------------------------------------------------------
%   Bz     -> Z component of the magnetic field measured in a regular planar grid
%   fs     -> Sampling frequency in 1/m
%   ----------------------------------------------------------------------------------------------


SHOWGRAPHS = 0; % Set this value to 1 to show frequency response plots, or to 0 otherwise.

[SIZEx, SIZEy] = size(Bz);
N1 = SIZEx;
N2 = SIZEy;

f1 = [0:N1 / 2, -(N1 / 2 - 1):-1] * fs / N1; %these freq. coordinates match the fft algorithm
f2 = [0:N2 / 2, -(N2 / 2 - 1):-1] * fs / N2;

ff1 = (-N1 / 2:N1 / 2 - 1) * fs / N1; %these freq. coordinates are more suitable to visualization
ff2 = (-N2 / 2:N2 / 2 - 1) * fs / N2;

[F2, F1] = meshgrid(f2+1e-30, f1+1e-30);
[FF2, FF1] = meshgrid(ff2+1e-30, ff1+1e-30);
kx = 2 * pi * F1;
ky = 2 * pi * F2;

kkx = 2 * pi * FF1;
kky = 2 * pi * FF2;

k = sqrt(kx.^2+ky.^2);
kk = sqrt(kkx.^2+kky.^2);

etx = -1i * kx ./ k; % calculate the filter frequency response associated with the x component
evx = -1i * kkx ./ kk;

ety = -1i * ky ./ k; % calculate the filter frequency response associated with the y component
evy = -1i * kky ./ kk;

e = fft2(Bz, N1, N2);

if SHOWGRAPHS
    figure
    surf(kkx, kky, imag(evx));
    figure
    surf(kkx, kky, imag(evy));
end

Bx = real(ifft2(e .* etx));
By = real(ifft2(e .* ety));
