function add_ROI(nROI, kwargs)
%add_ROI(nROI; 'ax', 'figure')
% adds a rectangle where the ROI is to a plot and labels it
% Parameters
% ----------
%     positional
%     ==========
%         ROI: cell
%         
%     keyword
%     =======
%         ax: Axes object (gca)
%         
arguments
    nROI cell
    kwargs.figure = gcf
    kwargs.ax = gca
end

fig = kwargs.figure;
ax = kwargs.ax;

set(fig, 'currentaxes', ax)

for i = 1 : size(nROI,2)
    iROI = nROI{i};
    [x, y, w, h] = get_mask_extent(iROI);
    rectangle('Position', [x y w h])
    t = text(double(x),double(y), ['#' num2str(i)]);
    set(t, 'Clipping', 'on')
end