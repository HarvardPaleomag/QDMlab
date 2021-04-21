function upward_continue(kwargs)

arguments
    kwargs.nFiles = 'none';
    kwargs.UC = 'none';
    kwargs.unit = 'T';
    kwargs.save {mustBeBoolean(kwargs.save)} = true;
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
            [b_uc] = UpCont(Btrunc, upDist, 1/step); %upward continue Bz map
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