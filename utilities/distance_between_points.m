function out = distance_between_points(data1, data2)
%[out] = distance_between_points(data1, data2; 'data2')
% calculates the distance between two points that get picked
arguments
    data1
    data2 = 'none'
end

nAx = 1;
if ~strcmp(data2, 'none')
    nAx = 2;
end

fig = figure('Name', 'QDM map', 'units', 'normalized', 'position', [0.1, 0.1, 0.4*nAx, 0.6]);
movegui(fig,'center')

ax1 = subplot(1,nAx,1);
QDM_figure(data1, 'ax', ax1)

if nAx == 2
    ax2 = subplot(1,nAx,2);
    QDM_figure(data2, 'ax', ax2)
end

ax1 = subplot(1,nAx,1);
[c1, r1, key] = ginput(1); %used to be c,l

for i = 1:nAx
    ax = subplot(1,nAx, i);
    hold on
    plot(ax, c1,r1, 'or')
end

ax2 = subplot(1,nAx,2);
[c2, r2, key] = ginput(1); %used to be c,l

for i = 1:nAx
    ax = subplot(1,nAx, i);
    hold on
    plot(ax, c2,r2, 'or')
end

c = [c1 c2];
r = [r1 r2];

out.p1 = [r(1), c(1)];
out.p2 = [r(2), c(2)];

out.distance = sqrt(diff(r)^2 + diff(c)^2);

axes = [];
for i = 1:nAx
    ax = subplot(1,nAx, i);
    hold on
    c = sort(c);
    r = sort(r);
    plot(ax, c,r, 'x-g');
    plot(ax, [c(1), c(1)],sort(r), 'x--g');
    plot(ax, c, [r(2), r(2)], 'x--g');

    text(mean(c)*1.05, mean(r), sprintf('d = %.1f px', out.distance))
    text(mean(c), r(2)*1.05, sprintf('dx = %.1f px', diff(c)), 'HorizontalAlignment', 'center')
    text(c(1)*0.95, mean(r), sprintf('dy = %.1f px', diff(r)), 'HorizontalAlignment', 'right')
    axes = [axes ax];
end

linkaxes(axes)
