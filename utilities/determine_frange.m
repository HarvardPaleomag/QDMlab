function [frange, findex] = determine_frange(folder)
arguments
    folder = 'none'
end

folder = automatic_input_ui__(folder, "single",1, "type",'dir');
[meanDataAll, freqAll] = load_mean_spectra(folder);

figure
hold on

for i = 1:size(meanDataAll, 1)
    plot(freqAll(i,:), meanDataAll(i,:), '.-', 'DisplayName', sprintf('run (%i)', i-1))
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
msg = sprintf('start1: %i    end1: %i', frange(1), frange(2));
logMsg('info',msg,1,0);
msg = sprintf('start2: %i    end2: %i', frange(3), frange(4));
logMsg('info',msg,1,0);
