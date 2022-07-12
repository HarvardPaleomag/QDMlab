function [idx] = pick_fcrop(meanData,freq)
%[idx] = pick_fcrop(meanData, freq)
    figure
    axis tight
    plot(freq, meanData, '.-');
    xlabel('contrast');
    ylabel('f (GHz');
    hold on
    idx = [0,0];
    for i = 1:2
        [x, ~] = ginput(1);
        [~, n] = min(abs(freq - x));
        idx(i) = n;
        plot(freq(n), meanData(n), 'xr')
    end
    plot(freq(idx(1):idx(2)), meanData(idx(1):idx(2)), 'ok-');
    msg = sprintf('cropping data between index %i : %i >> reduced datapoints to %i', idx(1), idx(2), diff(idx)+1);
    logMsg('debug',msg,1,0);
    msg = sprintf('New frequency Range: %i : %i', freq(idx(1)), freq(idx(2)));
    logMsg('info',msg,1,0);
    pause(2)
    close gcf
end