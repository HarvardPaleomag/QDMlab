function results = dipole_fit(kwargs)
%[results] = dipole_fit('filePath', 'fitOrder', 'xy', 'cropFactor', 'downSample', 'nRuns', 'quad', 'outputtrue', 'checkPlot', 'statsPlot', 'save', 'constrained', 'm0', 'hguess', 'minheight', 'maxheight', 'boxwidth', 'method', 'noise', 'SNR', 'AUTO', 'minTol', 'display', 'expData', 'dx', 'dy', 'imagefolder', 'sourceName')
%fitOrder is : 
%
%
% Parameters
% ----------
%     filePath: ['none']
%     fitOrder: [1]
%       highest order in the multipole expansion
%       1: dipole; 2: quadrupole; 8: octupole
%     xy: (int int) ['picker']
%       Location of the dipole in pixel
%       IF 'picker', lets you pick the location of the dipole
%     cropFactor: int [20]
%       Width of th box around 'xy' to be cropped
%     downSample: [1]
%       Option to downsample the data in order to speed up calculations
%     nRuns: [10]
%       Number of calculation runs
%     quad: []
%     save: [true]
%       saves data to file
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
%     STATISTSICS: int [0]
%     NSTAT: int [50]
%     display:
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
    kwargs.filePath = 'none';
    kwargs.fitOrder = 'none';
    kwargs.xy = 'picker';
    kwargs.cropFactor = 'none';
    kwargs.downSample = 1; % speedup
    kwargs.nRuns = 10;
    kwargs.quad = 1;
    kwargs.outputtrue {mustBeBoolean(kwargs.outputtrue)} = true;
    kwargs.checkPlot (1,1) {mustBeBoolean(kwargs.checkPlot)} = true;
    kwargs.statsPlot (1,1) {mustBeBoolean(kwargs.statsPlot)} = false;
    kwargs.save {mustBeBoolean(kwargs.save)} = 'none';
    
    kwargs.constrained {mustBeBoolean(kwargs.constrained)} = false; 
    kwargs.m0 = 1e-12;
    kwargs.hguess = 2.5e-5;
    kwargs.minheight = 0;
    kwargs.maxheight = 100e-6;
    kwargs.boxwidth = 100e-6;

    kwargs.method = 1; %0 = least squares, 1=Nelder-Mead
    kwargs.noise = 0; %0 = no noise
    kwargs.SNR = 0; %signal-to-noise ratio in dB
    kwargs.AUTO = 0; %automatically find dipole position from Bt map
    kwargs.minTol = 1;
    kwargs.display {mustBeBoolean(kwargs.display)} = false;
    
    kwargs.expData = 'none'; % loaded data passed -> no need to load data again
    kwargs.dx = false
    kwargs.dy = false
    
    kwargs.imagefolder = 'none';
    kwargs.sourceName = 'none'
    
end

% define defaults for the function
defaults = struct('fitOrder', 1, 'cropFactor', 20, 'save', true);

terms = [3,8,15];

% change random number generator
rng('shuffle')

set(0, 'DefaultFigureColormap', jet)

outFileName = 'DipoleInversions.txt';

% get the correct filePath - even if expData is passed
filePath = automatic_input_ui__(kwargs.filePath, 'type', 'file', 'single', true);

% get other parameters
kwargs = ask_arguments(kwargs, defaults);

if isa(kwargs.expData, 'struct')
    expData = kwargs.expData;
else
    msg = sprintf('loading data file:  %s', filePath);
    logMsg('info',msg,1,0);
%     fprintf('<>   loading data file:  %s\n', filePath);
    expData = load(filePath);
end


if strcmp(kwargs.xy, 'picker')
    [y,x] = pick_box2('expData', expData, 'point',true);
    XY = [x y];
else
    XY = kwargs.xy;
end

%% check NV distance (h)
if isfield(expData, 'h')
    h = expData.h;
else
    %ask user for the NV layer-sample distance
    msg = logMsg('input', 'NV-sample distance [5 µm]:', 0,0,'returnOnly', true);
    h = input(msg);
    if isempty(h)
        h=5e-6;
    end
end
% check if h/step is too large
if h > 0.1
    h = h * 1e-6;
end

%% check pixel size (step)
if isfield(expData, 'step')
    step = expData.step; % Pixel Size
else
    %ask user for the NV layer-sample distance
    msg = logMsg('input', 'Pixel Size [4.70 (µm)]:', 0,0,'returnOnly', true);
    step = input(msg);
    if isempty(step)
        step = 4.7e-6;
    end
end

if step > 1e-2
    step = step*1e-6;
end

%%
[bool, dataName,ledName] = is_B111(expData);
led = expData.(ledName);

[filePath, name, ext] = fileparts(filePath);

% downsample data
bData = downsample(downsample(expData.(dataName), kwargs.downSample)', kwargs.downSample)'; % Bz is assumed in T
bData = double(bData); % convert to double in case of single values

step = step * kwargs.downSample;

x = ((1:size(bData, 2)) - 1) * step;
y = ((1:size(bData, 1)) - 1) * step;
[X, Y] = meshgrid(x, y);

% in case the fitting is shown, reduce number of runs to 1
if kwargs.display
    kwargs.nRuns = 1;
end

j = round(XY(1)/kwargs.downSample);
i = round(XY(2)/kwargs.downSample);

%% CROP Data
% Transform lower/left index (i.e. from pick_sources_and_fit) to center location
% cropi, cropj are the cropping indices for x(i) and y(j) of the
% original map

% check if the optional parameters dx and dy are not 0
if all([kwargs.dx,kwargs.dy])
    % define cropj, cropi from dx, dy instead of CROPFACT
    xCrop = round([j j+kwargs.dx]);
    yCrop = round([i i+kwargs.dy]);
    i = i+kwargs.dy/2; j = j+kwargs.dx/2;
else
    xCrop = round(1+j+[-kwargs.cropFactor kwargs.cropFactor]);
    yCrop = round(1+i+[-kwargs.cropFactor kwargs.cropFactor]);
end

x0 = j * step;
y0 = i * step;

%adjust if the crop area falls outside the image
scansize = size(bData);

for p = 1:2
    if xCrop(p) < 1
        xCrop(p) = 1;
    end
    if yCrop(p) < 1
        yCrop(p) = 1;
    end
    if xCrop(p) > scansize(2)
        xCrop(p) = scansize(2);
    end
    if yCrop(p) > scansize(1)
        yCrop(p) = scansize(1);
    end
end

bDataCropped = bData(yCrop(1):yCrop(2), xCrop(1):xCrop(2));
Xc = X(yCrop(1):yCrop(2), xCrop(1):xCrop(2));
Yc = Y(yCrop(1):yCrop(2), xCrop(1):xCrop(2));
xc = x(xCrop(1):xCrop(2));
yc = y(yCrop(1):yCrop(2));

if kwargs.checkPlot
    figure
    imagesc(xc, yc, abs(bDataCropped));
    axis xy, axis equal, axis tight
    caxis([0, 1]*max(abs(caxis)));
    colormap(hot)
    colorbar
    title(sprintf('Cropped map'));
    drawnow
end

P00(1) = x0;
P00(2) = y0;
P00(3) = kwargs.hguess;

P = zeros(length(P00)+terms(kwargs.fitOrder), kwargs.nRuns);
fval = zeros(1, kwargs.nRuns);
fval2 = zeros(1, kwargs.nRuns);

for k = 1:kwargs.nRuns
    if kwargs.nRuns == 1
        P0 = P00; %+0.3*(rand(size(P00))-0.5).*P00;
    else
        P0 = P00 + 0.1 * (rand(size(P00)) - 0.5) .* P00;
    end

    options = optimset('TolX', 10^(floor(log10(kwargs.m0)) - 5), ...
        'TolFun', 10^(floor(log10(max(abs(bData(:))))) - 8), ...
        'MaxFunEvals', 6000, 'MaxIter', 2000, ...
        'Display', 'none');

    %% actual fitting
    if kwargs.method == 1
        
        if kwargs.constrained
            [P(1:3, k), fval2(k), exitflag, output] = fmincon(@(Pp) ...
                SourceFitMultiP8(Pp, Xc, Yc, bDataCropped, kwargs.display, ...
                kwargs.method, kwargs.quad, kwargs.fitOrder), P0, [], [], [], [], ...
                [x0 - kwargs.boxwidth, y0 - kwargs.boxwidth, kwargs.minheight], ...
                [x0 + kwargs.boxwidth, y0 + kwargs.boxwidth, kwargs.maxheight], [], options);
        else
            %% OLD unconstrained
            [P(1:3,k),fval2(k), exitflag, output] = fminsearch(@(Pp) ...
                SourceFitMultiP8(Pp,Xc,Yc,bDataCropped,kwargs.display,...
                kwargs.method, kwargs.quad, kwargs.fitOrder),P0,options);
        end

    else
        [P(1:3, k), fval2(k), resd, exitflag, output] = lsqnonlin(@(Pp) ...
            SourceFitMultiP8(Pp, Xc, Yc, bDataCropped, kwargs.display, ...
            kwargs.method, kwargs.quad, kwargs.fitOrder), P0, [], [], options);
    end

    %% calculate the model after fitting
    [resid, BzModel, M] = SourceFitMultiP8(P(1:3, k), Xc, Yc, bDataCropped, ...
                 kwargs.display, kwargs.method, kwargs.quad, kwargs.fitOrder);

    Mx = M(1);
    My = M(2);
    Mz = M(3);
    m = sqrt(Mx^2+My^2+Mz^2);

    theta = acosd(Mz/m);
    phi = atan2d(My, Mx);

    P(4, k) = m;
    %convert angles to inclination and declination
    P(5, k) = 90 - theta; % inc
    P(6, k) = phi - 90; % dec

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

    %Quadrupole moment
    if kwargs.fitOrder > 1
        P(7, k) = M(4);
        P(8, k) = M(5);
        P(9, k) = M(6);
        P(10, k) = M(7);
        P(11, k) = M(8);
    end

    %Octupole moment
    if kwargs.fitOrder > 2
        P(12, k) = M(9);
        P(13, k) = M(10);
        P(14, k) = M(11);
        P(15, k) = M(12);
        P(16, k) = M(13);
        P(17, k) = M(14);
        P(18, k) = M(15);
    end
   
    if size(P, 1) > 6
        Paux = [P(1:4, k)', 90 - P(5, k), P(6, k) + 90, P(7:terms(kwargs.fitOrder) + 3, k)'];
    else
        Paux = [P(1:4, k)', 90 - P(5, k), P(6, k) + 90];
    end

    [~, bModel] = SourceFitMultiP8(Paux, Xc, Yc, bDataCropped, 0, kwargs.method, kwargs.quad, kwargs.fitOrder);
    fval(k) = sqrt(sum(sum((bModel - bDataCropped).^2))/numel(bDataCropped)); %??????????????????????????
    
    %% display calculation and counter
    perc = round(k/kwargs.nRuns *50);
    strOut = [sprintf('<>                      FITTING: (%02i/%02i) [', k, kwargs.nRuns) repmat('*',1,perc) repmat(' ',1,50-perc),']'];
    % str to replace all previous characters
    strCR = repmat('\b',1,length(strOut));
    if k ~= 1; fprintf(strCR); end
    % display new line, replacing the old line
    fprintf(strOut);
    
    if k == kwargs.nRuns
        fprintf('\n')
    end
    
end

i0 = find(fval == min(fval));
fsort = sort(fval);

i = find(fval <= fsort(kwargs.minTol));
if numel(i) > 1
    msg = sprintf('averaging %d points\n', numel(i));
    logMsg('info',msg,1,0);
end
if numel(i) > 0.1 * numel(fval)
    msg = sprintf('Too many points are being averaged. Consider adjusting MINTOL parameter.');
    logMsg('warn',msg,1,0);
end

%% determine results from parameters
% moment
Popt = zeros(size(Paux));
Popt(4) = sum(P(4, i).*fval(i)) / sum(fval(i));
mopt = Popt(4);
% Inclination
Popt(5) = sum(P(5, i).*fval(i)) / sum(fval(i));
iopt = Popt(5);
% declination
Popt(6) = sum(P(6, i).*fval(i)) / sum(fval(i));
dopt = Popt(6); %x and y axis are reversed in SM - N is Y instead of x

% height
Popt(1) = sum(P(1, i).*fval(i)) / sum(fval(i));
Popt(2) = sum(P(2, i).*fval(i)) / sum(fval(i));
Popt(3) = sum(P(3, i).*fval(i)) / sum(fval(i));
hopt = Popt(3);

for kk = 7:terms(kwargs.fitOrder) + 3
    Popt(kk) = sum(P(kk, i).*fval(i)) / sum(fval(i));
end

% location
xopt = Popt(1);
yopt = Popt(2);
% fprintf('<>   RESULTS: \n');
msg = sprintf('M = %1.3d; I = %1.3f; D = %1.3f', mopt, iopt, dopt);
logMsg('result',msg,1,1);
msg = sprintf('h = %1.3d; x = %1.3d; y = %1.3d', hopt, xopt, yopt);
logMsg('result',msg,1,1);

%% calculate residuals
% parameters for model
Popt2 = [Popt(1:4), 90 - Popt(5), Popt(6) + 90, Popt(7:terms(kwargs.fitOrder) + 3)];
% calculate the model
[resid, bModel] = SourceFitMultiP8(Popt2, Xc, Yc, bDataCropped, 0, kwargs.method, kwargs.quad, kwargs.fitOrder);
% [residFull, bModelFull] = SourceFitMultiP8(Popt2, X, Y, bData, 0, kwargs.method, kwargs.quad, kwargs.fitOrder);

residuals = bModel - bDataCropped;
resids = sqrt(sum(sum(residuals.^2))/sum(sum(bDataCropped.^2)));

dipolarity = 1-(rms(residuals-mean(mean(residuals)))/rms(bDataCropped-mean(mean(bDataCropped))));
msg = sprintf('dipolarity parameter: %.2f', dipolarity);
logMsg('RESULT',msg,1,1);
    
nameext = [name, ext];

%% plotting / saving
if kwargs.checkPlot
    % whole map figure
    f = figure('Units', 'normalized', ...
               'Position',[0.1 0.2 0.8 0.4],'NumberTitle', 'off', 'Name', 'total map / LED');
           
    ax1 = subplot(1, 2, 1);
    QDM_figure(bData, 'ax', ax1);
    hold on
    rectPos = [xCrop(1), yCrop(1), xCrop(2)-xCrop(1), yCrop(2)-yCrop(1)];
    rectangle('Position', rectPos,...
        'EdgeColor', 'r', 'FaceColor', 'none', 'LineWidth', 0.7);
    
    ax2 = subplot(1, 2, 2);
    QDM_figure(led, 'ax', ax2, 'led',true, 'title', 'LED map');
    binning = detect_binning(expData);
    hold on
    rectPos = [xCrop(1)*binning, yCrop(1)*binning, (xCrop(2)-xCrop(1))*binning, (yCrop(2)-yCrop(1))*binning];
    rectangle('Position', rectPos,...
        'EdgeColor', 'r', 'FaceColor', 'none', 'LineWidth', 0.7);
    
    % data figure
    dataFig = figure();
    ax1 = subplot(2, 2, 1);
    imagesc(xc, yc, bDataCropped);
    axis xy, axis equal, axis tight;
    caxis([-1, 1]*max(abs(caxis)));
    colorbar
    title('Original Scan');
    ax2 = subplot(2, 2, 2);
    imagesc(xc, yc, bModel);
    axis xy, axis equal, axis tight;
    caxis([-1, 1]*max(abs(caxis)));
    colorbar
    title('Model Scan');
    ax3 = subplot(2, 2, 4);
    imagesc(xc, yc, residuals);
    axis xy, axis equal, axis tight;
    caxis([-1, 1]*max(abs(caxis)));
    colorbar
    title('Residuals');
    linkaxes([ax1, ax2, ax3]);
end


% create the outputs of the funtion
if hopt < 0
    dec = mod(180 - dopt, 360);
else
    dec = mod(360 - dopt, 360);
end

if kwargs.save
    % save figure
    saveas(dataFig, [filePath, '/Fit_', kwargs.sourceName, '_M', num2str(kwargs.fitOrder), '_x', num2str(round(XY(1))), 'y', num2str(round(XY(2))), '.png'])
    
    % add line to dipoleinversions.txt
    fid = fopen([filePath, '/', outFileName], 'r');
    header = (fid == -1);
    if ~header
        fclose(fid);
    end
    fid = fopen([filePath, '/', outFileName], 'a+t');
    if header
        fprintf(fid, 'File Name\tMoment\tInclination\tDeclination\tHeight\tDipolarity\r\n');
    end
    %Note a 180 rotation about y axis is imposed here
    fprintf(fid, '%s\t%1.5d\t%1.5d\t%1.5d\t%1.5d\t%1.5d\r\n', nameext, mopt, -iopt, dec, abs(hopt), dipolarity);
    fclose(fid);
end

if kwargs.statsPlot
    checkPlotFigure(P, fval, i, i0, mopt, iopt, dopt, hopt,xopt, yopt)
end

% create the outputs of the funtion
results = struct('dfile', filePath, 'm', mopt, 'inc', -iopt, 'dec', dec, ...
    'h', -hopt, 'res', resids, 'x',xopt,'y',yopt, 'residuals', residuals,...
    'data', bDataCropped, 'model', bModel, 'dipolarity', dipolarity, ...
    'xCrop', xCrop, 'yCrop', yCrop);
%     'Popt', Popt2, 'residFull',residFull, 'bModelFull',bModelFull % todo
%     save full model
end

function checkPlotFigure(P, fval, i, i0, mopt, iopt, dopt, hopt,xopt, yopt)
%checkPlotFigure(P, fval, i, i0, mopt, iopt, dopt, hopt, xopt, yopt)
    figure
    subplot(2,3,1)
    plot(P(1, :), fval, '.')
    title('Moment');
    hold on
    plot(P(1, i), fval(i), 'ko')
    plot(P(1, i0), fval(i0), 'go')
    plot([mopt, mopt], ylim, 'm--');
    hold off

    subplot(2,3,2)
    plot(P(2, :), fval, '.')
    title('Inclination');
    hold on
    plot(P(2, i), fval(i), 'ko')
    plot(P(2, i0), fval(i0), 'go')
    plot([iopt, iopt], ylim, 'm--');
    hold off

    subplot(2,3,3)
    plot(P(3, :), fval, '.')
    title('Declination');
    hold on
    plot(P(3, i), fval(i), 'ko')
    plot(P(3, i0), fval(i0), 'go')
    plot([dopt, dopt], ylim, 'm--');
    hold off

    subplot(2,3,4)
    plot(P(6, :), fval, '.')
    title('Height');
    hold on
    plot(P(6, i), fval(i), 'ko')
    plot(P(6, i0), fval(i0), 'go')
    plot([hopt, hopt], ylim, 'm--');
    hold off

    subplot(2,3,5)
    plot(P(4, :), fval, '.')
    title('X displacement');
    hold on
    plot(P(4, i), fval(i), 'ko')
    plot(P(4, i0), fval(i0), 'go')
    plot([xopt, xopt], ylim, 'm--');
    hold off

    subplot(2,3,6)
    plot(P(5, :), fval, '.')
    title('Y displacement');
    hold on
    plot(P(5, i), fval(i), 'ko')
    plot(P(5, i0), fval(i0), 'go')
    plot([yopt, yopt], ylim, 'm--');
    hold off
end

% 
% function results = dipole_fit(kwargs)
% %mOrder is : 
% %
% %
% % Parameters
% % ----------
% %     filePath: ['none']
% %     mOrder: [1]
% %       highest order in the multipole expansion
% %       1: dipole; 2: quadrupole; 8: octupole
% %     xy: (int int) ['picker']
% %       Location of the dipole in pixel
% %       IF 'picker', lets you pick the location of the dipole
% %     cropFactor: int [20]
% %       Width of th box around 'xy' to be cropped
% %     downSample: [1]
% %       Option to downsample the data in order to speed up calculations
% %     nRuns: [10]
% %       Number of calculation runs
% %     quad: []
% %     save: [true]
% %       saves data to file
% %     checkPlot: [false]
% % 
% %     incl: [-0]
% %     dec: [0]
% % 
% %     dx: int
% %         default: 0
% %         the width of the box that is cropped around the data
% %     dy: int
% %         default: 0
% %         the height of the box that is cropped around the data
% %     data: struct
% %         default: 0
% %         If you want to load data directly instead of reloading it from the
% %         file you can pass it after the 'data' keyword.
% %
% %     following parameters for costraining fit
% %     ========================================
% %     m0: double [1e-12]
% %     hguess: double [2.5e-5]
% %     minheight: double [0]
% %     maxheight: double [100e-6]
% %     boxwidth: double [100e-6]
% %
% %
% %     STATISTSICS: int [0]
% %     NSTAT: int [50]
% %     display:
% %
% % Returns: struct
% %     returns a structure with the following entries:
% %         dfile: file path
% %         m: fitted moment
% %         inc: fitted inclination
% %         dec: fitted declination
% %         h: fitted height (depth) of dipole
% %         res: residuals
% 
% arguments
%     kwargs.filePath = 'none';
%     kwargs.mOrder = 'none';
%     kwargs.xy = 'picker';
%     kwargs.cropFactor = 'none';
%     kwargs.downSample = 1; % speedup
%     kwargs.nRuns = 10;
%     kwargs.quad = 1;
%     kwargs.outputtrue {mustBeBoolean(kwargs.outputtrue)} = true;
%     kwargs.checkPlot {mustBeBoolean(kwargs.checkPlot)} = true;
%     kwargs.save {mustBeBoolean(kwargs.save)} = 'none';
% 
%     kwargs.m0 = 1e-12;
%     kwargs.hguess = 2.5e-5;
%     kwargs.minheight = 0;
%     kwargs.maxheight = 100e-6;
%     kwargs.boxwidth = 100e-6;
% 
%     kwargs.method = 1; %0 = least squares, 1=Nelder-Mead
%     kwargs.noise = 0; %0 = no noise
%     kwargs.SNR = 0; %signal-to-noise ratio in dB
%     kwargs.AUTO = 0; %automatically find dipole position from Bt map
%     kwargs.minTol = 1;
%     kwargs.statistics {mustBeBoolean(kwargs.statistics)} = false;
%     kwargs.nStats {mustBeInteger(kwargs.nStats)} = 50; % n iterations for statistics
%     kwargs.display {mustBeBoolean(kwargs.display)} = false;
%     
%     kwargs.expData = 'none'; % loaded data passed -> no need to load data again
%     kwargs.dx = false
%     kwargs.dy = false
%     
%     kwargs.imagefolder = 'none';
%     kwargs.sourceName = 'none'
%     
% end
% st = dbstack; 
% funcName = st.name; 
% 
% % define defaults for the function
% defaults = struct('mOrder', 1, 'cropFactor', 20, 'save', true);
% 
% terms = [3,8,15];
% 
% % change random number generator
% rng('shuffle')
% 
% set(0, 'DefaultFigureColormap', jet)
% 
% outFileName = 'DipoleInversions.txt';
% 
% % get the correct filePath - even if expData is passed
% filePath = automatic_input_ui__(kwargs.filePath, 'type', 'file', 'single', true);
% % get other parameters
% kwargs = ask_arguments(kwargs, defaults);
% 
% if isa(kwargs.expData, 'struct')
%     expData = kwargs.expData;
% else
%     msg = sprintf('loading data file:  %s', filePath);
%     logMsg('info',funcName,msg,1,0);
% %     fprintf('<>   loading data file:  %s\n', filePath);
%     expData = load(filePath);
% end
% 
% 
% if strcmp(kwargs.xy, 'picker')
%     [y,x] = pick_box2('expData', expData, 'point',true);
%     XY = [x y];
% else
%     XY = kwargs.xy;
% end
% 
% if exists_struct(expData, 'h')
%     h = expData.h;
% else
%     %ask user for the NV layer-sample distance
%     h = input('<>   INPUT << NV-sample distance [5 µm] >>: ');
%     if isempty(h)
%         h=5e-6;
%     end
% end
% 
% 
% step = expData.step; % Pixel Size
% 
% % check if h/step is too large
% if h > 0.1
%     h = h * 1e-6;
% end
% 
% if step > 1e-2
%     step = step*1e-6;
% end
% 
% [bool, dataName,ledName] = is_B111(expData);
% led = expData.(ledName);
% 
% [filePath, name, ext] = fileparts(filePath);
% 
% % downsample data
% bData = downsample(downsample(expData.(dataName), kwargs.downSample)', kwargs.downSample)'; % Bz is assumed in T
% bData = double(bData); % convert to double in case of single values
% 
% step = step * kwargs.downSample;
% 
% x = ((1:size(bData, 2)) - 1) * step;
% y = ((1:size(bData, 1)) - 1) * step;
% [X, Y] = meshgrid(x, y);
% 
% % in case the fitting is shown, reduce number of runs to 1
% if kwargs.display
%     kwargs.nRuns = 1;
% end
% 
% j = round(XY(1)/kwargs.downSample);
% i = round(XY(2)/kwargs.downSample);
% 
% %% CROP Data
% % Transform lower/left index (i.e. from pick_sources_and_fit) to center location
% % cropi, cropj are the cropping indices for x(i) and y(j) of the
% % original map
% 
% % check if the optional parameters dx and dy are not 0
% if all([kwargs.dx,kwargs.dy])
%     % define cropj, cropi from dx, dy instead of CROPFACT
%     cropj = round([j j+kwargs.dx]);
%     cropi = round([i i+kwargs.dy]);
%     i = i+kwargs.dy/2; j = j+kwargs.dx/2;
% else
%     cropj = round(1+j+[-kwargs.cropFactor kwargs.cropFactor]);
%     cropi = round(1+i+[-kwargs.cropFactor kwargs.cropFactor]);
% end
% 
% x0 = j * step;
% y0 = i * step;
% 
% %adjust if the crop area falls outside the image
% scansize = size(bData);
% 
% for p = 1:2
%     if cropj(p) < 1
%         cropj(p) = 1;
%     end
%     if cropi(p) < 1
%         cropi(p) = 1;
%     end
%     if cropj(p) > scansize(2)
%         cropj(p) = scansize(2);
%     end
%     if cropi(p) > scansize(1)
%         cropi(p) = scansize(1);
%     end
% end
% 
% bDataCropped = bData(cropi(1):cropi(2), cropj(1):cropj(2));
% Xc = X(cropi(1):cropi(2), cropj(1):cropj(2));
% Yc = Y(cropi(1):cropi(2), cropj(1):cropj(2));
% xc = x(cropj(1):cropj(2));
% yc = y(cropi(1):cropi(2));
% 
% if kwargs.checkPlot
%     figure
%     imagesc(xc, yc, abs(bDataCropped));
%     axis xy, axis equal, axis tight
%     caxis([0, 1]*max(abs(caxis)));
%     colormap(hot)
%     colorbar
%     title(sprintf('Cropped map'));
%     drawnow
% end
% 
% P00(1) = x0;
% P00(2) = y0;
% P00(3) = kwargs.hguess;
% 
% P = zeros(length(P00)+terms(kwargs.mOrder), kwargs.nRuns);
% fval = zeros(1, kwargs.nRuns);
% fval2 = zeros(1, kwargs.nRuns);
% 
% for k = 1:kwargs.nRuns
%     if kwargs.nRuns == 1
%         P0 = P00; %+0.3*(rand(size(P00))-0.5).*P00;
%     else
%         P0 = P00 + 0.1 * (rand(size(P00)) - 0.5) .* P00;
%     end
% 
%     options = optimset('TolX', 10^(floor(log10(kwargs.m0)) - 5), ...
%         'TolFun', 10^(floor(log10(max(abs(bData(:))))) - 8), ...
%         'MaxFunEvals', 6000, 'MaxIter', 2000, ...
%         'Display', 'none');
% 
%     %% actual fitting
%     if kwargs.method
%         [P(1:3, k), fval2(k), exitflag, output] = fmincon(@(Pp) ...
%             SourceFitMultiP8(Pp, Xc, Yc, bDataCropped, kwargs.display, ...
%             kwargs.method, kwargs.quad, kwargs.mOrder), P0, [], [], [], [], ...
%             [x0 - kwargs.boxwidth, y0 - kwargs.boxwidth, 2e-5], ...
%             [x0 + kwargs.boxwidth, y0 + kwargs.boxwidth, 3e-5], [], options);
%     else
%         [P(1:3, k), fval2(k), resd, exitflag, output] = lsqnonlin(@(Pp) ...
%             SourceFitMultiP8(Pp, Xc, Yc, bDataCropped, kwargs.display, ...
%             kwargs.method, kwargs.quad, kwargs.mOrder), P0, [], [], options);
%     end
% 
%     %% calculate the model after fitting
%     [resid, BzModel, M] = SourceFitMultiP8(P(1:3, k), Xc, Yc, bDataCropped, ...
%                  kwargs.display, kwargs.method, kwargs.quad, kwargs.mOrder);
% 
%     Mx = M(1);
%     My = M(2);
%     Mz = M(3);
%     m = sqrt(Mx^2+My^2+Mz^2);
% 
%     theta = acosd(Mz/m);
%     phi = atan2d(My, Mx);
% 
%     P(4, k) = m;
%     %convert angles to inclination and declination
%     P(5, k) = 90 - theta; % inc
%     P(6, k) = phi - 90; % dec
% 
%     %enforce range for parameters
%     change = 1;
%     while change
%         change = 0;
%         if P(5, k) > 90
%             P(5, k) = 180 - P(5, k); % i' = 180-i
%             P(6, k) = P(6, k) + 180; % d' = d+180
%             change = 1;
%         end
%         if P(5, k) < -90
%             P(5, k) = -180 - P(5, k); % i' = -180-i
%             P(6, k) = P(6, k) + 180; % d' = d+180
%             change = 1;
%         end
%         if P(6, k) < 0
%             P(6, k) = P(6, k) + 360; % d' = d+360
%             change = 1;
%         end
%         if P(6, k) >= 360
%             P(6, k) = P(6, k) - 360; % d' = d-360
%             change = 1;
%         end
%     end
% 
%     %Quadrupole moment
%     if kwargs.mOrder > 1
%         P(7, k) = M(4);
%         P(8, k) = M(5);
%         P(9, k) = M(6);
%         P(10, k) = M(7);
%         P(11, k) = M(8);
%     end
% 
%     %Octupole moment
%     if kwargs.mOrder > 2
%         P(12, k) = M(9);
%         P(13, k) = M(10);
%         P(14, k) = M(11);
%         P(15, k) = M(12);
%         P(16, k) = M(13);
%         P(17, k) = M(14);
%         P(18, k) = M(15);
%     end
%    
%     if size(P, 1) > 6
%         Paux = [P(1:4, k)', 90 - P(5, k), P(6, k) + 90, P(7:terms(kwargs.mOrder) + 3, k)'];
%     else
%         Paux = [P(1:4, k)', 90 - P(5, k), P(6, k) + 90];
%     end
% 
%     [~, bModel] = SourceFitMultiP8(Paux, Xc, Yc, bDataCropped, 0, kwargs.method, kwargs.quad, kwargs.mOrder);
%     fval(k) = sqrt(sum(sum((bModel - bDataCropped).^2))/numel(bDataCropped)); %??????????????????????????
%     
%     %% display calculation and counter
%     perc = round(k/kwargs.nRuns,2) *100;
%     strOut = [sprintf('<>   FITTING: (%02i/%02i) [', k, kwargs.nRuns) repmat('*',1,perc) repmat(' ',1,100-perc),']'];
%     % str to replace all previous characters
%     strCR = repmat('\b',1,length(strOut));
%     if k ~= 1; fprintf(strCR); end
%     % display new line, replacing the old line
%     fprintf(strOut);
%     
%     if k == kwargs.nRuns
%         fprintf('\n')
%     end
%     
% end
% 
% i0 = find(fval == min(fval));
% fsort = sort(fval);
% 
% i = find(fval <= fsort(kwargs.minTol));
% if numel(i) > 1
%     msg = sprintf('averaging %d points\n', numel(i));
%     logMsg('info',funcName,msg,1,0);
% end
% if numel(i) > 0.1 * numel(fval)
%     msg = sprintf('Too many points are being averaged. Consider adjusting MINTOL parameter.');
%     logMsg('warn',funcName,msg,1,0);
% end
% 
% %% determine results from parameters
% % moment
% Popt = zeros(size(Paux));
% Popt(4) = sum(P(4, i).*fval(i)) / sum(fval(i));
% mopt = Popt(4);
% % Inclination
% Popt(5) = sum(P(5, i).*fval(i)) / sum(fval(i));
% iopt = Popt(5);
% % declination
% Popt(6) = sum(P(6, i).*fval(i)) / sum(fval(i));
% dopt = Popt(6); %x and y axis are reversed in SM - N is Y instead of x
% 
% % height
% Popt(1) = sum(P(1, i).*fval(i)) / sum(fval(i));
% Popt(2) = sum(P(2, i).*fval(i)) / sum(fval(i));
% Popt(3) = sum(P(3, i).*fval(i)) / sum(fval(i));
% hopt = Popt(3);
% 
% for kk = 7:terms(kwargs.mOrder) + 3
%     Popt(kk) = sum(P(kk, i).*fval(i)) / sum(fval(i));
% end
% 
% % location
% xopt = Popt(1);
% yopt = Popt(2);
% fprintf('<>   RESULTS: \n');
% fprintf('<>      M = %1.3d (min = %1.3d); I = %1.3f (min = %1.3f); D = %1.3f (min = %1.3f) \n', mopt, P(4, i0), iopt, P(5, i0), dopt, P(6, i0));
% fprintf('<>      h = %1.3d (min = %1.3d); x = %1.3d (min = %1.3d); y = %1.3d (min = %1.3d)\n', hopt, P(3, i0), xopt, P(1, i0), yopt, P(2, i0));
% 
% 
% %% calculate residuals
% % parameters for model
% Popt2 = [Popt(1:4), 90 - Popt(5), Popt(6) + 90, Popt(7:terms(kwargs.mOrder) + 3)];
% % calculate the model
% [resid, bModel] = SourceFitMultiP8(Popt2, Xc, Yc, bDataCropped, 0, kwargs.method, kwargs.quad, kwargs.mOrder);
% residex = bModel - bDataCropped;
% resids = sqrt(sum(sum(residex.^2))/sum(sum(bDataCropped.^2)));
% 
% nameext = [name, ext];
% 
% %% plotting / saving
% if kwargs.checkPlot
%     % whole map figure
%     f = figure('Units', 'normalized', ...
%                'Position',[0.1 0.2 0.8 0.4],'NumberTitle', 'off', 'Name', 'total map / LED');
%            
%     ax1 = subplot(1, 2, 1);
%     QDM_figure(bData, 'ax', ax1);
%     hold on
%     rectPos = [cropj(1), cropi(1), cropj(2)-cropj(1), cropi(2)-cropi(1)];
%     rectangle('Position', rectPos,...
%         'EdgeColor', 'r', 'FaceColor', 'none', 'LineWidth', 0.7);
%     
%     ax2 = subplot(1, 2, 2);
%     QDM_figure(led, 'ax', ax2, 'led',true, 'title', 'LED map');
%     binning = detect_binning(expData);
%     hold on
%     rectPos = [cropj(1)*binning, cropi(1)*binning, (cropj(2)-cropj(1))*binning, (cropi(2)-cropi(1))*binning];
%     rectangle('Position', rectPos,...
%         'EdgeColor', 'r', 'FaceColor', 'none', 'LineWidth', 0.7);
%     
%     % data figure
%     dataFig = figure();
%     ax1 = subplot(2, 2, 1);
%     imagesc(xc, yc, bDataCropped);
%     axis xy, axis equal, axis tight;
%     caxis([-1, 1]*max(abs(caxis)));
%     colorbar
%     title('Original Scan');
%     ax2 = subplot(2, 2, 2);
%     imagesc(xc, yc, bModel);
%     axis xy, axis equal, axis tight;
%     caxis([-1, 1]*max(abs(caxis)));
%     colorbar
%     title('Model Scan');
%     ax3 = subplot(2, 2, 4);
%     imagesc(xc, yc, residex);
%     axis xy, axis equal, axis tight;
%     caxis([-1, 1]*max(abs(caxis)));
%     colorbar
%     title('Residuals');
%     linkaxes([ax1, ax2, ax3]);
% end
% 
% 
% % create the outputs of the funtion
% if hopt < 0
%     dec = mod(180 - dopt, 360);
% else
%     dec = mod(360 - dopt, 360);
% end
% 
% if kwargs.save
%     % save figure
%     saveas(dataFig, [filePath, '/Fit_', name, '_M', num2str(kwargs.mOrder), '_x', num2str(round(XY(1))), 'y', num2str(round(XY(2))), '.png'])
%     
%     % add line to dipoleinversions.txt
%     fid = fopen([filePath, '/', outFileName], 'r');
%     header = (fid == -1);
%     if ~header
%         fclose(fid);
%     end
%     fid = fopen([filePath, '/', outFileName], 'a+t');
%     if header
%         fprintf(fid, 'File Name\tMoment\tInclination\tDeclination\tHeight\tResiduals\r\n');
%     end
%     %Note a 180 rotation about y axis is imposed here
%     fprintf(fid, '%s\t%1.5d\t%1.5d\t%1.5d\t%1.5d\t%1.5d\r\n', nameext, mopt, -iopt, dec, abs(hopt), resids);
%     fclose(fid);
% end
% 
% % create the outputs of the funtion
% results = struct('dfile', filePath, 'm', mopt, 'inc', -iopt, 'dec', dec, 'h', -hopt, 'res', resids, 'x',xopt,'y',yopt,'residuals',residex);
% end
% 
% function checkPlotFigure(P, fval, i, i0, mopt, iopt, dopt, hopt,xopt, yopt)
% 
% figure
% subplot(2,3,1)
% plot(P(1, :), fval, '.')
% title('Moment');
% hold on
% plot(P(1, i), fval(i), 'ko')
% plot(P(1, i0), fval(i0), 'go')
% plot([mopt, mopt], ylim, 'm--');
% hold off
% 
% subplot(2,3,2)
% plot(P(2, :), fval, '.')
% title('Inclination');
% hold on
% plot(P(2, i), fval(i), 'ko')
% plot(P(2, i0), fval(i0), 'go')
% plot([iopt, iopt], ylim, 'm--');
% hold off
% 
% subplot(2,3,3)
% plot(P(3, :), fval, '.')
% title('Declination');
% hold on
% plot(P(3, i), fval(i), 'ko')
% plot(P(3, i0), fval(i0), 'go')
% plot([dopt, dopt], ylim, 'm--');
% hold off
% 
% subplot(2,3,4)
% plot(P(6, :), fval, '.')
% title('Height');
% hold on
% plot(P(6, i), fval(i), 'ko')
% plot(P(6, i0), fval(i0), 'go')
% plot([hopt, hopt], ylim, 'm--');
% hold off
% 
% subplot(2,3,5)
% plot(P(4, :), fval, '.')
% title('X displacement');
% hold on
% plot(P(4, i), fval(i), 'ko')
% plot(P(4, i0), fval(i0), 'go')
% plot([xopt, xopt], ylim, 'm--');
% hold off
% 
% subplot(2,3,6)
% plot(P(5, :), fval, '.')
% title('Y displacement');
% hold on
% plot(P(5, i), fval(i), 'ko')
% plot(P(5, i0), fval(i0), 'go')
% plot([yopt, yopt], ylim, 'm--');
% hold off
% end
% 
