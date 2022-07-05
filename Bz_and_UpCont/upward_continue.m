function ucMaps = upward_continue(kwargs)
%[ucMaps] = upward_continue('nFiles', 'UC', 'unit', 'save')
%
% Parameters
% ----------
%   nFiles: ('none')
%   UC: ('none')
%   unit: ('T')
%   save: (true)
% 
% Returns
% ----------
%   ucMaps:

arguments
    kwargs.nFiles = 'none';
    kwargs.UC = 'none';
    kwargs.unit = 'T';
    kwargs.save {mustBeBoolean(kwargs.save)} = true;
end

nFiles = automatic_input_ui__(kwargs.nFiles, 'type', 'file', ...
    'title', 'Pick a magnetic field map file');

if strcmp(kwargs.UC, 'none')
    msg = sprintf('Enter UC distances (micron) with [ ] around them: ');
    msg = logMsg('INPUT',msg,0,0,'returnOnly',true);
    kwargs.UC = input(msg);
end

UC = kwargs.UC * 1e-6; %upward continuations in m

[filePath, name, ~] = fileparts(nFiles{:});

msg = sprintf('loading: %s', nFiles{:});
logMsg('info',msg,1,0);

expData = load(nFiles{:});

% if a B111dataToPlot.mat file is passed, the step, h variables is not
% saved -> need to define it.
if exist('step', 'var') == 0
    step = 4.68e-6; % pixel size in (m)
    h = 5e-6; % NV layer thickness (m) ?
end

[bool, dataName,ledName] = is_B111(expData);

Bdata = expData.(dataName);
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

        fig = QDM_figure(Bdata, ...
            'cbTitle', sprintf('B_z (%s)', kwargs.unit), ...
            'title', sprintf('QDM data UC = %.1f micron', upDist * 1e6),...
            'axis', 'off');

        if kwargs.save
            fileName = strrep(name, '_uc0.mat', '');
            fileName = [fileName sprintf('_uc%.0f', upDist * 1e6)];
            fileName = check_suffix(fileName);
            
%             fileName = [name, '_uc', num2str(upDist * 1e6), '.mat'];
            
            dataOut = expData;
            dataOut.([dataName '_uc0']) = Btrunc;
            dataOut.(dataName) = b_uc;
            dataOut.h = h;
            dataOut.step = step;
            dataOut.Bt = Bt;
            dataOut.fileName = fileName;
            dataOut.filePath = filePath;
            ucMaps{end+1} = dataOut;
            
            msg = sprintf('saving: %s and %s', fileName, strrep(fileName,'mat','png'));
            logMsg('info',msg,1,0);
            
            saveas(fig, fullfile(filePath, strrep(fileName,'mat','png')));
            save(fullfile(filePath, fileName), '-struct', 'dataOut', '-mat'); % why no 'Bt',
            

        end
    end
end

end

function Btrunc = makeEvenArray(Bdata)
%[Btrunc] = makeEvenArray(Bdata)
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