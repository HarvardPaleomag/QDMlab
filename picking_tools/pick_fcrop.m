function [idx] = pick_fcrop(meanData,freq)
    figure
    axis tight
    plot(freq, meanData, '.-');
    xlabel('contrast');
    ylabel('f (GHz');
    hold on
    idx = [0,0];
    for i = 1:2
        [x, ~] = ginput(1);
        abs(meanData - x)
        [~, n] = min(abs(freq - x));
        idx(i) = n;
        plot(freq(n), meanData(n), 'xr')
    end
    plot(freq(idx(1):idx(2)), meanData(idx(1):idx(2)), 'ok-')
    pause(5)
    close gcf
end