function [nROI, coordinates] = pick_box(data, kwargs)
%[nROI, coordinates] = pick_box(data; 'led', 'closeFig', 'returnCoordinates', 'title', 'n')
% 
% positional parameters
%     data: QDM/LED data
% 
% optional positional parameters
%     n: int
%     number of rectangles to pick
% 
% optional parameters:
%     led: bool
%         if true: colorscale 'bone' will be used for plotting
%     close_fig: bool
%         if true: figure will be closed after pick
% 
% Returns: cell
%     a cell with n entries. Each consisting of [x1,x2,y1,y2], where (x1,x2)

arguments
    data
    kwargs.led (1,1) {mustBeBoolean(kwargs.led)} = false
    kwargs.closeFig (1,1) {mustBeBoolean(kwargs.closeFig)} = false
    kwargs.returnCoordinates (1,1) {mustBeBoolean(kwargs.returnCoordinates)} = false
    kwargs.title = 'pick boxes'
    kwargs.n = 'none'
    kwargs.std = 2
end

% data = filter_hot_pixels(data);


if kwargs.led == 1
    fig = QDM_figure(data, 'led', true, 'title', [kwargs.title,' (ESC to exit)']);
else
    fig = QDM_figure(data, 'title', [kwargs.title,' (ESC to exit)'], 'std',kwargs.std);
end

figTitle = 'Pick Sources (ESC to exit)';

title(figTitle)
nROI = {};
coordinates = {};
n = 1;

while n
    iROI = drawrectangle(gca, 'Label', num2str(n), 'LabelVisible', 'on', ...
                         'HandleVisibility', 'off', 'InteractionsAllowed', 'none');
    
    if iROI.Position

        loc = iROI.Position;

        x0 = int16(loc(1));
        dx = int16(loc(3));
        x1 = x0 + dx;

        y0 = int16(loc(2));
        dy = int16(loc(4));
        y1 = y0 + dy;

        % create box for plotting
        box = [[x0, x0, x1, x1, x0]; [y0, y1, y1, y0, y0]];
        hold on
        plot(box(1, :), box(2, :), '-', 'linewidth', 2);
        drawnow


        %% coordinates
        % save all points you continue getting
        % rounded and negative values -> 0
        x0 = round(x0); y0 = round(y0); dx = round(dx); dy = round(dy);
        msg = sprintf('creating coordinates of box #%i lower left = (%i,%i) dx:%i, dy:%i', n, x0, y0, dx, dy);
        logMsg('info',msg,1,0);
        coordinates{end+1} = max(round([x0, y0, dx, dy]), 0);
        
        % ROI
        iMask = zeros(size(data));
        iMask(y0:y1,x0:x1)=1;
        m = limit_mask(iMask);
        msg = sprintf('creating mask for box #%i (%ix%i : %i pixel)', n, size(m,2), size(m,1), numel(m));
        logMsg('info',msg,1,0);
        nROI{end+1} = iMask;
        
        if n == kwargs.n
            n = false;
            continue
        end
        
        n = n + 1;
    else
        n = false;
        continue
    end
end

if kwargs.closeFig
    close(fig)
end