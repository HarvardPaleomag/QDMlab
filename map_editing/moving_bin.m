function binnedData = moving_bin(data, win, kwargs)
% bins data using a moving window. The window type can be selected.
%
% Parameters
% ----------
%     required
%     ========
%     data
%     win
%     
%     keyword
%     =======
%     winType
%     

arguments
    data double
    win double
    kwargs.winType {mustBeMember(kwargs.winType, ['boxcar', 'sine', 'tria', 'hamming', 'han'])} = 'boxcar'
end

tStart = tic;
[w,  wSum] = window(win, 'winType', kwargs.winType);

binnedData = zeros(size(data));
sizeX = size(data,2); sizeY = size(data,1);

for r = win+1:sizeY-win
    for c = win+1:sizeX-win
        dWin = data(r-win:r+win, c-win:c+win, :);        
        meanD = sum(dWin.*w, [1 2]);
        meanD = meanD / wSum;
        meanD = squeeze(meanD);
        binnedData(r,c,:) = meanD;
    end
end
fprintf('<>      INFO: finished moving window binning: %.1f s\n', toc(tStart)');
end

%% Window FUNCTIONS
function [w, wSum] = window(win, kwargs)
% creates a matrix with weights according to the win size and winType.
    arguments
        win
        kwargs.winType = 'boxcar'
    end
    
    N = 2*win+1;
    w = ones(N,N);
    w0 = [];    

    if strcmp(kwargs.winType, 'boxcar')
        wSum = sum(w, 'all');
        return
    elseif strcmp(kwargs.winType, 'tria')
        for n = 1:2*N
            if n > 2*win+1
                break
            end
            wn = 1 - abs((n-(N+1)/2) / ((N+1)/2));
            w0 = [w0 wn];
        end
    elseif strcmp(kwargs.winType, 'sine')
        for n = 1:2*N+1
            if n > 2*win+1
                break
            end
            wn = sin( pi * n/(N+1));
            w0 = [w0 wn];
        end
    elseif strcmp(kwargs.winType, 'hamming') | strcmp(kwargs.winType, 'hann')
        if strcmp(kwargs.winType, 'hamming') 
            a0 = 0.53836;
        elseif strcmp(kwargs.winType, 'hann')
            a0 = 0.5;
        end
        
        for n = 1:2*N
            if n > 2*win+1
                break
            end
            
            wn = a0 - (1-a0) * cos( 2*pi * n/(N+1));
            w0 = [w0 wn];
        end
    end

    w = w .* w0;
    w = w .* w0';
    wSum = sum(w, 'all');
end