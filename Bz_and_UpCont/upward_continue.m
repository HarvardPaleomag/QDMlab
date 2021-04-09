function upward_continue(kwargs)

arguments
    kwargs.nFiles = 'none';
    kwargs.UC = 'none';
    kwargs.unit = 'T';
    kwargs.rotate = 0;
    kwargs.medFilter = 0;
    kwargs.save {mustBeBoolean(kwargs.save)} = true;
    kwargs.cal = 1; %field calibration factor
    kwargs.alpha = 0; %rotation of the diamond lattice axes around z-axis
    kwargs.beta = 0; %rotation of the image axes around z-axis
end

nFiles = automatic_input_ui__(kwargs.nFiles, 'type', 'file', ...
    'title', 'Pick a magnetic field map file');

if strcmp(kwargs.UC, 'none')
    kwargs.UC = input('<>   Enter UC distances (micron) with [ ] around them: ');
end

UC = kwargs.UC * 1e-6; %upward continuations in m

[filepath, name, ~] = fileparts(nFiles{:});

expData = load(nFiles{:});

% if a B111dataToPlot.mat file is passed, the step, h variables is not
% saved -> need to define it.
if exist('step', 'var') == 0
    step = 4.68e-6; % pixel size in (m)
    h = 5e-6; % NV layer thickness (m) ?
end

if is_B111(expData)
    Bdata = expData.B111ferro;
else
    Bdata = expData.Bz;
end

Btrunc = makeEvenArray(Bdata);
ucMaps = {};

if size(UC, 2) ~= 0
    for n = 1:size(UC, 2)
        upDist = UC(n);

        if upDist == 0
            [b_uc] = Btrunc;
        else
            [b_uc] = calculate_UC(Btrunc, upDist, 1/step); %upward continue Bz map
        end

        [by, bx] = MITBxByFromBz(b_uc, 1/step);

        %show figures
        Bt = sqrt(bx.^2+by.^2+b_uc.^2);
        Bdata = b_uc; %Conserve kwargs.units
        h = h + upDist;

        fig = QDM_figure(Bdata, 'cbTitle', sprintf('B_z (%s)', kwargs.unit), 'axis', 'off');

        if kwargs.save
            fileName = [name, '_uc', num2str(upDist * 1e6), '.mat'];

            saveas(fig, [filepath, '/', name, '_uc', num2str(upDist * 1e6), '.png'])
            if exists_struct(expData, 'newLED')
                newLED = expData.newLED;
                corners = expData.corers;

                if is_B111(expData)
                    B111ferro = Bdata;
                    save(fullfile(filepath, fileName), 'B111ferro', 'Bt', 'h', 'step', 'newLED', 'corners', '-mat');
                else
                    Bz = Bdata;
                    save(fullfile(filepath, fileName), 'Bz', 'Bt', 'h', 'step', 'newLED', 'corners', '-mat');
                end
            else
                if is_B111(expData)
                    B111ferro = Bdata;
                    save(fullfile(filepath, fileName), 'B111ferro', 'h', 'step', '-mat'); % why no 'Bt',
                else
                    Bz = Bdata;
                    save(fullfile(filepath, fileName), 'Bz', 'h', 'step', '-mat'); % why no 'Bt',
                end
            end
            disp(fprintf(['<>   INFO: Saved: ', name, '_uc', num2str(upDist * 1e6), '.mat...\n']))
        end
    end
end

end

function Btrunc = makeEvenArray(Bdata)
%make even sized arrays
if mod(size(Bdata, 1), 2)
    if mod(size(Bdata, 2), 2)
        Btrunc = Bdata(2:end, 2:end);
    else
        Btrunc = Bdata(2:end, 1:end);
    end
else
    if mod(size(Bdata, 2), 2)
        Btrunc = Bdata(1:end, 2:end);
    else
        Btrunc = Bdata(1:end, 1:end);
    end
end
end

function [Bout] = calculate_UC(data, dh, fs)
% Parameters
% ----------
%       data: double
%           component map
%       dh: double
%           distance to upward continue the field map
%       fs: double
%           Sampling frequency in 1/m

EXPAND = 1;
% Algorithm refines sampling in the frequency domain by
% this factor (i.e., zero padding). Increase value (should be
% an even number) if the number of points in Bz map is small.

FOVERSAMPL = 2;

if EXPAND
    Bino = data;

    data = [zeros(size(data)), zeros(size(data)), zeros(size(data)); ...
        zeros(size(data)), data, zeros(size(data)); ...
        zeros(size(data)), zeros(size(data)), zeros(size(data))];
end

[sizeX, sizeY] = size(data');
Nx = sizeX * FOVERSAMPL;
Ny = sizeY * FOVERSAMPL;

% Adjust sampling frequency accordingly
% fs= fs * FOVERSAMPL;

% These freq. coordinates match MATLAB's fft algorithm
f2 = [0:Ny / 2, -(Ny / 2 - 1):-1] * fs / Ny;
f1 = [0:Nx / 2, -(Nx / 2 - 1):-1] * fs / Nx;

% Generate spatial frequency grid. A tiny 'epsilon' is added to
% each freq. variable so as to avoid division by zero warning.
[F2, F1] = meshgrid(f2, f1);

kx = 2 * pi * F1;
ky = 2 * pi * F2;
k = sqrt(kx.^2+ky.^2);

% Calculate the filter frequency response associated with the x component
et = exp(-dh*k);

% Compute FFT of the field map
e = fft2(data', Nx, Ny);

% Calculate x component
Bout = ifft2(e.*et, Nx, Ny, 'symmetric')';

% Crop matrices to get rid of zero padding
Bout = Bout(1:Ny/FOVERSAMPL, 1:Nx/FOVERSAMPL);

if EXPAND
    Bout = Bout(size(Bino, 1)+1:2*size(Bino, 1), size(Bino, 2)+1:2*size(Bino, 2));
end
end