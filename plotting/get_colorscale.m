function [cLim] = get_colorscale(data, method, kwargs, filter)
%[cLim] = get_colorscale(data; 'method', 'symmetric', 'mustBe', 'std', 'nOutlier')

arguments
    data
    method = 'std'
    kwargs.symmetric = 1
    kwargs.mustBe = 'none'
    filter.std = 0
    filter.nOutlier = 0
end

% Set the remaining axes properties

d = reshape(data,[numel(data), 1]);
d = d(~isnan(d));
d = sort(d);
d = d(1+filter.nOutlier:end-filter.nOutlier);

cLim = [min(d), max(d)];
switch method
    case 'std'
        med = median(abs(d), 'all', 'omitnan');
        st = std(d, [], 'all', 'omitnan');
        mx = max(d, [], 'all', 'omitnan');
        mn = min(d, [], 'all', 'omitnan');

        if (med + filter.std * st) > max(abs([mx,mn]))
            msg = sprintf('Clim values exceeds min/max');
            logMsg('debug',msg,1,0);
        elseif ~all(data(~isnan(data)) > 0, 'all')
            msg = sprintf('setting Clim: +-%.3e, according to: median (%.3e) + %i*std (%.3e)', med+filter.std*st, med,filter.std, st);
            logMsg('debug',msg,1,0);
            cLim = [-1, 1]*(med + filter.std * st);
        else
            msg = sprintf('setting Clim: %.3e:%.3e, according to: median (%.3e) +- %i*std (%.3e)', ...
                         med-filter.std*st, med+filter.std*st, med,filter.std, st);
            logMsg('info',msg,1,0);
            CLim = [med - filter.std * st, med + filter.std * st];
        end
    case 'fit'
        [p,S] = polyfit(1:numel(d),d,3);
        warning('off','MATLAB:polyfit:RepeatedPointsOrRescale')

        x = 1 : numel(d);
        [yFit,delta] = polyval(p,x,S);
        
        
        if kwargs.symmetric
            cLim = [-1 1] * (max(abs(yFit)) + filter.std * delta(1));
        else
            cLim = [min(yFit), max(yFit)] + [-1 1] * filter.std * delta(1);
        end
        
        msg = sprintf('setting Clim: %.2e:%.2e, according to: polyfit(3) ± %i std', ...
                     cLim(1), cLim(2),filter.std);
        logMsg('info',msg,1,0);
end

if strcmp(kwargs.mustBe,'pos')
    if cLim(1) < 0
        cLim(1) = 0;
    end
elseif strcmp(kwargs.mustBe,'neg')
    if cLim(2) > 0
        cLim(2) = 0;
    end
end

end