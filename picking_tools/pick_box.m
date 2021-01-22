function nRects = pick_box(data, varargin)
% pick_box lets you pick the ceter of an image and a box around it
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

p = inputParser;
str_or_char = @(x) isstring(x) | ischar(x);

addRequired(p, 'data');
addParameter(p, 'led', 'false', @islogical);
addParameter(p, 'closeFig', false, @islogical);
addParameter(p, 'returnCoordinates', false, @islogical);
parse(p, data, varargin{:});

led = p.Results.led;
closeFig = p.Results.closeFig;
retCoord = p.Results.returnCoordinates;

if led == 1
%     fig = LEDimage(data);
    fig = figure;
    imagesc(data);
else
    fig = QDM_figure(data);
end
figTitle = 'Pick Sources (ESC to exit)';

title(figTitle)
nRects = {};
n = 1;

while n
    %    rect = getrect(gca());
    iRect = drawrectangle(gca, 'Label', num2str(n), 'LabelVisible', 'on', ...
                         'HandleVisibility', 'off', 'InteractionsAllowed', 'none');

    if iRect.Position

        loc = iRect.Position;

        x0 = loc(1);
        dx = loc(3);
        x1 = x0 + dx;

        y0 = loc(2);
        dy = loc(4);
        y1 = y0 + dy;

        % create box for plotting
        box = [[x0, x0, x1, x1, x0]; [y0, y1, y1, y0, y0]];
        hold on
        plot(box(1, :), box(2, :), '-', 'linewidth', 2);
        drawnow

        % save all points you continue getting
        % rounded and negative values -> 0
        if retCoord
            fprintf('<>      creating coordinates of box #%i (%i,%i) dx:%i, dy:%i', n, x0, y0, dx, dy)
            nRects{end+1} = max(round([x0, y0, dx, dy]), 0);
        else
            iMask = createMask(iRect);
            m = limit_mask(iMask);
            fprintf('<>      creating mask for box #%i (%ix%i : %i pixel)\n', n, size(m,2), size(m,1), numel(m))
            nRects{end+1} = iMask;
        end

        %
        n = n + 1;
    else
        n = false;
        continue
    end
end

if closeFig
    close(fig)
end