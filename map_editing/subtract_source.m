function expData = subtract_source(kwargs)
% This script takes an input Bz map, asks for a box, crops to that box, and
% outputs Bz and Bt maps, along with the accessory parameters
%
% Parameters
% ----------
%     filePath: str ['none']
%     save: bool [true]
%     fitOrder: 1, 2 or 3 [1]
%       1=dipole; 2=quadrupole; 3=octapole
%
% Returns
% -------
%     Map: double
%       after subtraction

arguments
    kwargs.filePath = 'none'
    kwargs.save = false;
    kwargs.fitOrder = 1;
    kwargs.checkPlot = true;
end

filePath = automatic_input_ui__(kwargs.filePath, 'type', 'file', 'title', 'Pick a magnetic field map file');
[filePath, fileName, ext] = fileparts(filePath);

expData = load(filePath{:});
expData.filePath = filePath;

[b111, dataName, ~] = is_B111(expData);

if b111
    error('B111 file selected, please use Bz file')
end

bData = expData.(dataName);
expData.([dataName, '_original']) = bData;

[expData, row, col] = crop_map('filePath', expData, 'save', false, 'checkPlot', false);

residualMap = FitMoment(kwargs.fitOrder, filePath{:}, [col(1, 1), row(1, 1)], [col(2, 1), row(2, 1)], 2, 0);

ic = 1;
for i = col(1):col(2)
    jc = 1;
    for j = row(1):row(2)
        bData(j, i) = -residualMap(jc, ic);
        jc = jc + 1;
    end
    ic = ic + 1;
end

expData.(dataName) = bData;

if kwargs.checkPlot
    QDM_figure(bData);
end

% calculate new Btotal
[by, bx] = MITBxByFromBz(bData, 1/expData.step);
expData.Bt = sqrt(by.^2+bx.^2+bData.^2);

if kwargs.save
    save(fullfile(filepath, [fileName(1:end - 1), num2str(str2num(fileName(end)) + 1), '.mat']), '-struct', 'expData');
end

end