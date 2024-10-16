function demag_behavior_plot(results, kwargs)
%demag_behavior_plot(results; 'steps', 'stepUnit', 'dataUnit', 'parameter', 'ylabel', 'led', 'mean')
%coercivity_result_plot(results; 'steps', 'stepUnit', 'led', 'mean')
% plots results from estimate_coercivity
% 
% Parameters
% ----------
%     positional
%     ==========
%         results: struct
%             results structure from 'estmate_coercivity'
%     keyword
%     =======
%         steps: bool (0)
%         led: bool (0)
%
arguments
    results struct
    kwargs.steps  (1,:) double = false
    kwargs.stepUnit  (1,:) = 'mT'
    kwargs.dataUnit  (1,:) = 'T'
    kwargs.parameter = 'pPixels'
    kwargs.ylabel = 'norm. n(+)pixel'
    kwargs.led  (1,1) {mustBeMember(kwargs.led, [1, 0])} = 0
    kwargs.mean (1,1) {mustBeBoolean(kwargs.mean)} = false
end
    
[nMasks, nFiles] = size(results.nFiles);

if isequal(kwargs.steps, false)
    steps = 1:nFiles;
else
    steps = kwargs.steps;
    if iscell(steps)
        steps = cell2mat(steps);
    end
end

if size(steps) ~= size(results.(kwargs.parameter), 2)
    s = size(steps,2);
    s2 = size(results.(kwargs.parameter),2);
    error('WARNING: number of steps (%i) does not match number of results (%i)', s, s2);
end

rows = fix(nFiles/3)+ double(mod(nFiles,3)>0);

figure('units','normalized', 'outerposition',[0,0,3*(0.8/rows),0.8]);

axes = [];
mx = 0;
mn = 0;
means = rand(1, nFiles);
stds = rand(1, nFiles);

for j = 1:nFiles
    % load the data of this file
    iFileData = results.transDatas{j};
    data(abs(iFileData) > mean(iFileData, 'all','omitnan') + 100 * std(iFileData, 0, 'all','omitnan') ) = nan;
    iFileData = iFileData - median(iFileData, 'all','omitnan');

    %% get max and min data values -> global max min
    if max(iFileData, [], 'all') > mx
        mx = max(iFileData, [], 'all', 'omitnan');
    end
    if min(iFileData, [], 'all') < mn
        mn = min(iFileData, [], 'all', 'omitnan');
    end
    
    means(j) = mean(iFileData, 'all', 'omitnan');
    stds(j) = std(iFileData, 0, 'all', 'omitnan');
    
    fileName = results.nFiles{1, j};
    fNameSplit = split(fileName,filesep);
    
    title = fNameSplit{end-2};
    
    if ~ all(steps == 1:nFiles,'all')
        title = [sprintf('%.1f %s', steps(j), kwargs.stepUnit)];
    end

    ax = subplot(rows, 3, j);
    QDM_figure(iFileData, 'ax', ax, 'title', title, 'unit', kwargs.dataUnit);
    hold(ax, 'on');
    axes = [axes ax];
    
    for i = 1:nMasks
        iMask = results.nMasks{i};
        nROI = results.nROI{i};
        visboundaries(iMask, 'lineWidth', 0.7);
        
        [x, y, w, h] = get_mask_extent(nROI);
        rectangle(ax, 'Position', [x y w h])
        t = text(ax, double(x),double(y), ['#' num2str(i)]);
        set(t, 'Clipping', 'on')
    end
end
linkaxes(axes)

for ax = axes
    lim = mean(means)+5*mean(stds);
    set(ax,'CLim',[-1 1] * lim);
end

figure
ay = subplot(1, 1, [1]);
hold(ay, 'on');

for i = 1:nMasks
    errorbar(ay, steps, results.(kwargs.parameter)(i, :, 1) / results.(kwargs.parameter)(i, 1, 1), ...
             results.(kwargs.parameter)(i, :, 2) /results.(kwargs.parameter)(i, 1, 1), ...
             'o-', 'DisplayName', num2str(i))
end

if kwargs.mean
    mn = mean(results.(kwargs.parameter)(:, :, 1) ./ results.(kwargs.parameter)(:, 1, 1)); 
    st = std(results.(kwargs.parameter)(:, :, 1) ./ results.(kwargs.parameter)(:, 1, 1)); 
    errorbar(ax, steps, mn(:,:,1), st(:,:,1), 'k.--', 'DisplayName', 'mean', 'LineWidth', 1)
end
ylabel(ay, kwargs.ylabel)
legend(ay)
h = plot(ay, [min(steps), max(steps)], [0.5, 0.5], '--', 'color', '#7F7F7F', 'DisplayName', '');
% the following line skip the name of the previous plot from the legend
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
%% additional LED plot
if kwargs.led
    axes = [];
    figure 
    for j = 1:nFiles
    % load the data of this file
        iFileLed = results.transLeds{j};

        ax = subplot(rows-1, 3, j);
        imagesc(ax, iFileLed)
        colormap('gray');
        axis equal, axis tight, axis xy

        hold on
        title(ax, step);
        axes = [axes ax];
    
        for i = 1:nMasks
            iMask = results.nMasks{i};
            nROI = results.nROI{i};
            
            iMask = re_bin(iMask, iFileLed);
            nROI = re_bin(nROI, iFileLed);
            
            visboundaries(iMask, 'lineWidth', 0.7);

            [x, y, w, h] = get_mask_extent(nROI);
            rectangle(ax, 'Position', [x y w h])
            t = text(ax, double(x),double(y), ['#' num2str(i)]);
            set(t, 'Clipping', 'on')
        end
    end
linkaxes(axes)
end
