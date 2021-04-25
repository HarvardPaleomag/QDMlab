function results = dipole_fit(kwargs)
%kwargs.mOrder is highest order in the multipole expansion: 2 is quadrupole; 8 is
%octupole
%
%
% Parameters
% ----------
%     dataFile: ['none']
%     mOrder: [1]
%     xy: [[0, 0]]
%     cropFactor: []
%     downSample: [1] %WHY?
%     nRuns: []
%     quad: []
%     outputtrue: [true]
%     checkPlot: [false]
% 
%     incl: [-0]
%     dec: [0]
% 
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
%     m0: double [1e-12]
%     hguess: double [2.5e-5]
%     minheight: double [0]
%     maxheight: double [100e-6]
%     boxwidth: double [100e-6]
%
%
%     STATISTSICS:
%         default: 0
%     NSTAT:
%         default: 50
%     kwargs.display:
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

arguments
    kwargs.dataFile = 'none';
    kwargs.mOrder = 'none';
    kwargs.xy = [0, 0];
    kwargs.cropFactor = 100;
    kwargs.downSample = 1; % speedup
    kwargs.nRuns = 10;
    kwargs.quad = 1;
    kwargs.outputtrue {mustBeBoolean(kwargs.outputtrue)} = true;
    kwargs.checkPlot {mustBeBoolean(kwargs.checkPlot)} = false;

    kwargs.incl = -0;
    kwargs.dec = 0;

    kwargs.m0 = 1e-12;
    kwargs.hguess = 2.5e-5;
    kwargs.minheight = 0;
    kwargs.maxheight = 100e-6;
    kwargs.boxwidth = 100e-6;

    kwargs.method = 1; %0 = least squares, 1=Nelder-Mead
    kwargs.noise = 0; %0 = no noise
    kwargs.SNR = 0; %signal-to-noise ratio in dB
    kwargs.AUTO = 0; %automatically find dipole position from Bt map
    kwargs.TERMS = [3, 8, 15];
    kwargs.MINTOL = 1;
    kwargs.statistics {mustBeBoolean(kwargs.statistics)} = false;
    kwargs.nStats {mustBeInteger(kwargs.nStats)} = 50; % n iterations for statistics
    kwargs.display {mustBeBoolean(kwargs.display)} = false;

end

defaults = kwargs;
kwargs = ask_arguments(kwargs, defaults);

theta0 = 90 - kwargs.incl;
phi0 = kwargs.dec + 90; %SM y axis is the one pointing northy
XY = kwargs.xy;
terms = [3,8,15];

% change random number generator
rng('shuffle')

set(0, 'DefaultFigureColormap', jet)

outFileName = 'DipoleInversions.txt';

if kwargs.statistics
    counter = kwargs.nStats;
else
    counter = 1;
end

dataFile = automatic_input_ui__(kwargs.dataFile, 'type', 'file', 'single', true);

expData = load(dataFile);
[filePath, name, ext] = fileparts(dataFile);

[~, dataName, ~] = is_B111(expData);

while counter %get rid of statistics

    if kwargs.statistics
        fprintf('<>   Iteration %02i/%02i ...\n', kwargs.nStats-counter+1, kwargs.nStats);
    end

    step = 0;
    
    % downsample data
    Bdata = downsample(downsample(expData.(dataName), kwargs.downSample)', kwargs.downSample)'; % Bz is assumed in T
    Bdata = double(Bdata); % convert to double in case of single values
    
    step = step * kwargs.downSample;

    x = ((1:size(Bdata, 2)) - 1) * step;
    y = ((1:size(Bdata, 1)) - 1) * step;
    [X, Y] = meshgrid(x, y);

    if ~kwargs.statistics
        statsFigure = figure();
    end

    if kwargs.noise
        noise = randn(size(Bdata)) * 0.5 * std(std(Bdata));
        noise = sqrt(10^(-SNR / 10)*var(Bdata(:))) * randn(size(Bdata));
        Bdata = Bdata + noise;
    end

    if kwargs.display
        kwargs.nRuns = 1;
    end

    j = round(XY(1)/kwargs.downSample);
    i = round(XY(2)/kwargs.downSample);

    x0 = j * step;
    y0 = i * step;
    cropj = round(1+j+[-kwargs.cropFactor, kwargs.cropFactor]);
    cropi = round(1+i+[-kwargs.cropFactor, kwargs.cropFactor]);

    %adjust if the crop area falls outside the image
    scansize = size(Bdata);

    for p = 1:2
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

    BdataCropped = Bdata(cropi(1):cropi(2), cropj(1):cropj(2));
    Xc = X(cropi(1):cropi(2), cropj(1):cropj(2));
    Yc = Y(cropi(1):cropi(2), cropj(1):cropj(2));
    xc = x(cropj(1):cropj(2));
    yc = y(cropi(1):cropi(2));

    imagesc(xc, yc, abs(BdataCropped));
    axis xy, axis equal, axis tight
    caxis([0, 1]*max(abs(caxis)));
    colormap(hot)
    colorbar
    title(sprintf('Cropped map'));
    drawnow

    P00(1) = x0;
    P00(2) = y0;
    P00(3) = kwargs.hguess;

    drawnow
    P = zeros(length(P00)+terms(kwargs.mOrder), kwargs.nRuns);
    fval = zeros(1, kwargs.nRuns);
    fval2 = zeros(1, kwargs.nRuns);

    for k = 1:kwargs.nRuns
        if kwargs.nRuns == 1
            P0 = P00; %+0.3*(rand(size(P00))-0.5).*P00;
        else
            %P0=P00+0.1*(rand(size(P00))-0.5).*[P00(1) 20 40 P00(4:6)];
            P0 = P00 + 0.1 * (rand(size(P00)) - 0.5) .* P00;
        end

        options = optimset('TolX', 10^(floor(log10(kwargs.m0)) - 5), ...
            'TolFun', 10^(floor(log10(max(abs(Bdata(:))))) - 8), ...
            'MaxFunEvals', 6000, 'MaxIter', 2000, ...
            'Display', 'none');
        
        %% actual fitting
        if kwargs.method
            [P(1:3, k), fval2(k), exitflag, output] = fmincon(@(Pp) ...
                SourceFitMultiP8(Pp, Xc, Yc, BdataCropped, kwargs.display, ...
                kwargs.method, kwargs.quad, kwargs.mOrder), P0, [], [], [], [], ...
                [x0 - kwargs.boxwidth, y0 - kwargs.boxwidth, 2e-5], ...
                [x0 + kwargs.boxwidth, y0 + kwargs.boxwidth, 3e-5], [], options);
        else
            [P(1:3, k), fval2(k), resd, exitflag, output] = lsqnonlin(@(Pp) SourceFitMultiP8(Pp, Xc, Yc, BdataCropped, kwargs.display, kwargs.method, kwargs.quad, kwargs.mOrder), P0, [], [], options);
        end
        
        %% calculate the model after fitting
        [resid, BzModel, M] = SourceFitMultiP8(P(1:3, k), Xc, Yc, BdataCropped, kwargs.display, kwargs.method, kwargs.quad, kwargs.mOrder);
        
        Mx = M(1);
        My = M(2);
        Mz = M(3);
        m = sqrt(Mx^2+My^2+Mz^2);

        theta = acosd(Mz/m);
        phi = atan2d(My, Mx);

        P(4, k) = m;
        %convert angles to inclination and declination
        P(5, k) = 90 - theta;
        P(6, k) = phi - 90;

        %enforce range for parameters
        change = 1;
        while change
            change = 0;
            if P(5, k) > 90
                P(5, k) = 180 - P(5, k); % i' = 180-i
                P(6, k) = P(6, k) + 180; % d' = d+180
                change = 1;
            end
            if P(5, k) < -90
                P(5, k) = -180 - P(5, k); % i' = -180-i
                P(6, k) = P(6, k) + 180; % d' = d+180
                change = 1;
            end
            if P(6, k) < 0
                P(6, k) = P(6, k) + 360; % d' = d+360
                change = 1;
            end
            if P(6, k) >= 360
                P(6, k) = P(6, k) - 360; % d' = d-360
                change = 1;
            end
        end

        if kwargs.mOrder > 1
            %Quadrupole moment
            P(7, k) = M(4);
            P(8, k) = M(5);
            P(9, k) = M(6);
            P(10, k) = M(7);
            P(11, k) = M(8);
        end

        if kwargs.mOrder > 2
            %Octupole moment
            P(12, k) = M(9);
            P(13, k) = M(10);
            P(14, k) = M(11);
            P(15, k) = M(12);
            P(16, k) = M(13);
            P(17, k) = M(14);
            P(18, k) = M(15);
        end

        if kwargs.display
            fprintf('Moment: %.4e\n', P(4))
            fprintf('Inclination: %.1f\n', P(5))
            fprintf('Declination: %.1f\n', P(6))
            fprintf('Height: %.4e\n', P(3))
        end

        fprintf('..%0d..(%0d)  ', k, output.iterations);

        if ~mod(k, 10)
            fprintf('\n');
        end

        if size(P, 1) > 6
            Paux = [P(1:4, k)', 90 - P(5, k), P(6, k) + 90, P(7:terms(kwargs.mOrder) + 3, k)'];
        else
            Paux = [P(1:4, k)', 90 - P(5, k), P(6, k) + 90];
        end

        [resid, Bzmodel] = SourceFitMultiP8(Paux, Xc, Yc, BdataCropped, 0, kwargs.method, kwargs.quad, kwargs.mOrder);
        fval(k) = sqrt(sum(sum((Bzmodel - BdataCropped).^2))/numel(BdataCropped));
    end

    i0 = find(fval == min(fval));
    fsort = sort(fval);

    %i=find(fval<=min(fval)+MINTOL*(max(fval)-min(fval)));
    i = find(fval <= fsort(MINTOL));
    fprintf('\n--- Averaging %d points ---\n', numel(i))
    if numel(i) > 0.1 * numel(fval)
        disp('Too many points are being averaged. Consider adjusting MINTOL parameter.')
    end


    Popt = zeros(size(Paux));
    Popt(4) = sum(P(4, i).*fval(i)) / sum(fval(i));
    mopt = Popt(4); %*1000
    fprintf('(min = %1.3d)\n', P(4, i0));

    Popt(5) = sum(P(5, i).*fval(i)) / sum(fval(i));
    iopt = Popt(5);
    fprintf('(min = %1.3f)\n', P(5, i0));

    Popt(6) = sum(P(6, i).*fval(i)) / sum(fval(i));
    dopt = Popt(6); %x and y axis are reversed in SM - N is Y instead of x
    mod(360-dopt, 360)
    fprintf('(min = %1.3f)\n', P(6, i0));

    Popt(1) = sum(P(1, i).*fval(i)) / sum(fval(i));
    Popt(2) = sum(P(2, i).*fval(i)) / sum(fval(i));
    Popt(3) = sum(P(3, i).*fval(i)) / sum(fval(i));
    hopt = Popt(3);
    fprintf('(min = %1.3d)\n', P(3, i0));

    for kk = 7:terms(kwargs.mOrder) + 3
        Popt(kk) = sum(P(kk, i).*fval(i)) / sum(fval(i));
    end

    xopt = Popt(1);
    fprintf('(min = %1.3d)\n', P(1, i0));

    yopt = Popt(2);
    fprintf('(min = %1.3d)\n', P(2, i0));

    Popt2 = [Popt(1:4), 90 - Popt(5), Popt(6) + 90, Popt(7:terms(kwargs.mOrder) + 3)];
    [resid, Bzmodel] = SourceFitMultiP8(Popt2, Xc, Yc, BdataCropped, 0, kwargs.method, kwargs.quad, kwargs.mOrder);
    residex = Bzmodel - BdataCropped;

    if ~kwargs.statistics
        figure
        subplot(2, 2, 1);
        imagesc(xc, yc, BdataCropped);
        axis xy, axis equal, axis tight;
        caxis([-1, 1]*max(abs(caxis)));
        colorbar
        title('Original Scan');
        subplot(2, 2, 2);
        imagesc(xc, yc, Bzmodel);
        axis xy, axis equal, axis tight;
        caxis([-1, 1]*max(abs(caxis)));
        colorbar
        title('Model Scan');
        subplot(2, 2, 4);
        imagesc(xc, yc, residex);
        axis xy, axis equal, axis tight;
        caxis([-1, 1]*max(abs(caxis)));
        colorbar
        title('Residuals');
    end

    saveas(gcf, [filePath, '/Fit_', name, '_M', num2str(kwargs.mOrder), '_x', num2str(round(XY(1))), 'y', num2str(round(XY(2))), '.png'])

    %resids=sqrt(sum(sum(residex.^2))/numel(residex))
    resids = sqrt(sum(sum(residex.^2))/sum(sum(BdataCropped.^2)))

    nameext = [name, ext];
    if kwargs.statistics
        outputtrue = 1;
    end

    if outputtrue
        fid = fopen([filePath, '/', outFileName], 'r');
        header = (fid == -1);
        if ~header
            fclose(fid);
        end
        fid = fopen([filePath, '/', outFileName], 'a+t');
        if header
            fprintf(fid, 'File Name\tMoment\tInclination\tDeclination\tHeight\tResiduals\r\n');
        end
        %Note a 180 rotation about y axis is imposed here
        if hopt < 0
            fprintf(fid, '%s\t%1.5d\t%1.5d\t%1.5d\t%1.5d\t%1.5d\r\n', nameext, mopt, -iopt, mod(180 - dopt, 360), -hopt, resids);
        else
            fprintf(fid, '%s\t%1.5d\t%1.5d\t%1.5d\t%1.5d\t%1.5d\r\n', nameext, mopt, -iopt, mod(360 - dopt, 360), hopt, resids);
        end
        fclose(fid);
    end

    counter = counter - 1;

end

% create the outputs of the funtion
if hopt < 0
    results = struct('dfile', INFILE, 'm', mopt, 'inc', -iopt, 'dec', mod(180 - dopt, 360), 'h', -hopt, 'res', resids);
else
    results = struct('dfile', INFILE, 'm', mopt, 'inc', -iopt, 'dec', mod(360 - dopt, 360), 'h', -hopt, 'res', resids);
end
end

function checkPlotFigure(P, fval)

figure
plot(P(1, :), fval, '.')
title('Moment');
hold on
plot(P(1, i), fval(i), 'r.')
plot(P(1, i0), fval(i0), 'c.')
plot([mopt, mopt], ylim, 'm--');
hold off

figure
plot(P(2, :), fval, '.')
title('Inclination');
hold on
plot(P(2, i), fval(i), 'r.')
plot(P(2, i0), fval(i0), 'c.')
plot([iopt, iopt], ylim, 'm--');
hold off

figure
plot(P(3, :), fval, '.')
title('Declination');
hold on
plot(P(3, i), fval(i), 'r.')
plot(P(3, i0), fval(i0), 'c.')
plot([dopt, dopt], ylim, 'm--');
hold off

figure
plot(P(6, :), fval, '.')
title('Height');
hold on
plot(P(6, i), fval(i), 'r.')
plot(P(6, i0), fval(i0), 'c.')
plot([hopt, hopt], ylim, 'm--');
hold off

figure
plot(P(4, :), fval, '.')
title('X displacement');
hold on
plot(P(4, i), fval(i), 'r.')
plot(P(4, i0), fval(i0), 'c.')
plot([xopt, xopt], ylim, 'm--');
hold off

figure
plot(P(5, :), fval, '.')
title('Y displacement');
hold on
plot(P(5, i), fval(i), 'r.')
plot(P(5, i0), fval(i0), 'c.')
plot([yopt, yopt], ylim, 'm--');
hold off
end
