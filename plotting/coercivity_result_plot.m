function coercivity_result_plot(results, kwargs)
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
    kwargs.steps  (1,:) = 0
    kwargs.led  (1,1) {mustBeMember(kwargs.led, [1, 0])} = 0

end
    
[nMasks, nFiles] = size(results.nFiles);

if iscell(kwargs.steps)
    steps = kwargs.steps;
else
    steps = 1:nFiles;
    
end
if size(steps) ~= size(results.pPixels, 2)
    s = size(steps,2);
    s2 = size(results.pPixels,2);
    disp(['<>   WARNING: number of steps (' num2str(s) ') does not match number of results (' num2str(s2) ')' ])
    return
end

rows = fix(nFiles/3)+ double(mod(nFiles,3)>0)+1;

if iscell(steps)
    steps = cell2mat(steps);
end

figure
axes = [];
mx = 0;
mn = 0;

for j = 1:nFiles
    % load the data of this file
    iFileData = results.transDatas{j};
    data(abs(iFileData) > nanmean(iFileData, 'all') + 100 * nanstd(iFileData, 0, 'all') ) = nan;
    iFileData = iFileData - nanmedian(iFileData, 'all');

    %% get max and min data values -> global max min
    if max(iFileData, [], 'all') > mx
        mx = nanmax(iFileData, [], 'all');
    end
    if min(iFileData, [], 'all') < mn
        mn = nanmin(iFileData, [], 'all');
    end

    fileName = results.nFiles{1, j};
    fNameSplit = split(fileName,filesep);
    step = fNameSplit(end-2);

    ax = subplot(rows, 3, j);
    imagesc(ax, iFileData)
    axis equal, axis tight, axis xy

    hold on
    title(ax, step);
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
    set(ax,'CLim',[-1 1] * mx/2);
end

% figure
for i = 1:nMasks
    ax = subplot(rows, 3, [rows*3-2 rows*3-1 rows*3]);
    hold on
    errorbar(steps, results.pPixels(i, :, 1) /results.pPixels(i, 1, 1), ...
        results.pPixels(i, :, 2) /results.pPixels(i, 1, 1), ...
        'o-', 'DisplayName', num2str(i))
    ylabel('norm. n(+)pixel')
    legend
end
set(gcf,'position',[250,250,3*200,rows*150])


%% additional LED plot
if kwargs.led
    axes = [];
    figure 
    for j = 1:nFiles
    % load the data of this file
        iFileLed = results.transLeds{j};
    % 
    %     %% get max and min data values -> global max min
    %     if max(iFileData, [], 'all') > mx
    %         mx = nanmax(iFileData, [], 'all');
    %     end
    %     if min(iFileData, [], 'all') < mn
    %         mn = nanmin(iFileData, [], 'all');
    %     end
    % 
    %     fileName = results.nFiles{1, j};
    %     fNameSplit = split(fileName,filesep);
    %     step = fNameSplit(end-2);

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