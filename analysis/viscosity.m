function viscData = viscosity(kwargs)

arguments
    kwargs.visc = 'none';
    kwargs.nonVisc = 'none';
    kwargs.checkPlot (1,1) {mustBeBoolean(kwargs.checkPlot)} = false

end

viscousFile = automatic_input_ui__(kwargs.visc, 'type', 'file', 'single', true);
viscousFile = check_suffix(viscousFile);

nonViscousFile = automatic_input_ui__(kwargs.nonVisc, 'type', 'file', 'single', true);
nonViscousFile = check_suffix(nonViscousFile);

%% align images
[viscousFile_, name, ext] = fileparts(viscousFile);
[nonViscousFile_, name, ext] = fileparts(nonViscousFile);

fprintf('<>   loading and transforming data')
[transformedData, nFiles] = get_transformed_maps({nonViscousFile_, viscousFile_}, ...
    'fileName', name, 'checkPlot', kwargs.checkPlot);

%% load data
viscousData = transformedData(viscousFile).targetData;
nonViscousData = transformedData(nonViscousFile).refData;
refData = load(nonViscousFile);
%%
delta = viscousData - nonViscousData;

if kwargs.checkPlot
    viscFigure = figure('Units', 'normalized', ...
               'Position',[0.05 0.2 0.9 0.25], 'Name', 'viscosity');
    ax1 = subplot(1,3,1);
    QDM_figure(viscousData,'ax', ax1, 'title', 'viscous');
    ax2 = subplot(1,3,2);
    QDM_figure(nonViscousData,'ax', ax2, 'title', 'non viscous');
    ax3 = subplot(1,3,3);
    QDM_figure(delta,'ax', ax3, 'title', 'delta');
    linkaxes([ax1 ax2 ax3]);
end

viscData.viscousData = viscousData;
viscData.nonViscousData = nonViscousData;
viscData.viscosityMap = delta;

[~, dataName, ledName] = is_B111(refData)

viscData.LED = refData.(ledName);

if exists_struct(refData, 'laser')
    viscData.laser = refData.laser;
else
    viscData.laser = zeros(size(refData.(ledName),2));
end

fName = [name '_viscosity', '.mat'];
fprintf('<>   SAVING data to %s', fName);
save(fullfile(viscousFile_, fName), '-struct', 'viscData')