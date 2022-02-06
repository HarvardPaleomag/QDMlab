function [frange, findex] = determine_frange(folder)
arguments
    folder = 'none'
end

folder = automatic_input_ui__(folder, "single",1, "type",'dir');

runFiles = dir(fullfile(folder, 'run_*.mat'));

for i = 1:size(runFiles)
    runFiles(i).data = load(fullfile(runFiles(i).folder, ... 
        runFiles(i).name), 'disp1', 'disp2', 'freqList');
end

figure
hold on

meanDataAll = [];
freqAll = [];
for i = 1:size(runFiles)
    plot(runFiles(i).data.freqList, ...
        [runFiles(i).data.disp1, runFiles(i).data.disp2], '.:', ...
        'displayName', runFiles(i).name)
    meanDataAll = [meanDataAll; [runFiles(i).data.disp1, runFiles(i).data.disp2]];
    freqAll = [freqAll; runFiles(i).data.freqList];
end
legend()

meanData = mean(meanDataAll);
freq = mean(freqAll);

frange = [];
findex = [];
for i = 1:4
    [x, ~] = ginput(1);
    [~, n] = min(abs(freq - x));
    frange = [frange freq(n)];
    findex = [findex n];
    plot(freq(n), meanData(n), 'o', 'DisplayName', sprintf('P(%i)', i),...
        'LineWidth',3)
    if ~mod(i,2)
        plot(freq(findex(i-1):findex(i)),...
            meanData(findex(i-1):findex(i)), 'k-', ...
            'LineWidth',2,...
            'HandleVisibility','off')

    end
end

