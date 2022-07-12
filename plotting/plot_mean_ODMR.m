function plot_mean_ODMR(folder, kwargs)
%plot_mean_ODMR(folder; 'closeAll', 'polarity', 'xlim')
% 
% Parameters
% ----------
%   folder:
%   closeAll: (0)
%   polarity: (1)
%   xlim: (false)
% 
% Returns
% ----------

arguments
    folder
    kwargs.closeAll = 0
    kwargs.polarity = 1
    kwargs.xlim = false
end
    if kwargs.closeAll
        close all
    end
    
    idx = kwargs.polarity;
    [d,f] = load_mean_spectra(folder);
    l = plot(f(idx,:)/1e9, d(idx,:)/max(d(idx,:)), '.-');
    
    if ~kwargs.xlim
        xlim([min(f(idx,:)), max(f(idx,:))]/1e9)
    end
end
