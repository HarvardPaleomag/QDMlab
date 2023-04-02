function [Bout] = UpCont(Bin, dh, fs)
% Upward continue a component map to a higher elevation.
%
% INPUTS:
% - Bin: A 2D array representing the component map to be upward continued.
%   The array should have dimensions (Ny, Nx), where Ny is the number of
%   rows and Nx is the number of columns.
%
% - dh: The distance to upward continue the field map, in the same units as
%   the sampling interval. For example, if the sampling interval is meters,
%   then dh should be specified in meters as well.
%
% - fs: The sampling frequency of the component map, in units of 1/m. For
%   example, if the component map is sampled every 1 meter, then fs should
%   be 1.
%
% OUTPUTS:
% - Bout: A 2D array representing the upward continued component map. The
%   array has the same dimensions as the input array Bin.
%
% NOTES:
% - This function uses the Fourier transform to apply a frequency-domain
%   filter to the input component map. The filter is designed to account
%   for the effect of elevation on the field values.
% - The function assumes that the input component map is sampled uniformly
%   in space, and that the elevation is constant across the entire map.
% - The function uses zero-padding and oversampling in the Fourier domain
%   to improve the accuracy of the upward continuation. The amount of
%   oversampling can be adjusted by changing the value of FOVERSAMPL.
% - The function can be slow for large input maps or high values of dh.



% Set algorithm parameters
EXPAND = 1;     % Whether to expand the input map by mirroring the edges
FOVERSAMPL = 2; % The amount of oversampling to use in the Fourier domain

% Expand input map if necessary
if EXPAND
    Bino = Bin; % an even number) if the number of points in Bz map is small.

    Bin = [zeros(size(Bin)), zeros(size(Bin)), zeros(size(Bin)); ...
        zeros(size(Bin)), Bin, zeros(size(Bin)); ...
        zeros(size(Bin)), zeros(size(Bin)), zeros(size(Bin))];
end

% Compute size and sampling frequency of input map
[SIZEx, SIZEy] = size(Bin');
Nx = SIZEx * FOVERSAMPL;
Ny = SIZEy * FOVERSAMPL;
%fs= fs * FOVERSAMPL;                        % Adjust sampling frequency accordingly

% Compute frequency coordinates
f2 = [0:Ny / 2, -(Ny / 2 - 1):-1] * fs / Ny; % These freq. coordinates match MATLAB's fft algorithm
f1 = [0:Nx / 2, -(Nx / 2 - 1):-1] * fs / Nx;

[F2, F1] = meshgrid(f2, f1); % Generate spatial frequency grid. A tiny 'epsilon' is added to
% each freq. variable so as to avoid division by zero warning.
kx = 2 * pi * F1;
ky = 2 * pi * F2;
k = sqrt(kx.^2+ky.^2);

% Compute filter frequency response
et = exp(-dh*k); % Calculate the filter frequency response associated with the x component

% Compute Fourier transform of input map
e = fft2(Bin', Nx, Ny); % Compute FFT of the field map

% Apply filter in frequency domain and compute inverse Fourier transform
Bout = ifft2(e.*et, Nx, Ny, 'symmetric')'; % Calculate x component

% Crop output map to remove zero-padding
Bout = Bout(1:Ny/FOVERSAMPL, 1:Nx/FOVERSAMPL); % Crop matrices to get rid of zero padding

if EXPAND
    Bout = Bout(size(Bino, 1)+1:2*size(Bino, 1), size(Bino, 2)+1:2*size(Bino, 2));
end
