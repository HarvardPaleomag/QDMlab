function nMasks = pick_area(data)
%[nMasks] = pick_area(data)
% Function lets you pick several areas of an image and returns N masks, one
% for each area. The Mask is an array with the same dimensiona as the
% original image.
% 
% returns:
%     cell of n logical arrays one for each mask

figure
imagesc(data)
axis xy; axis equal; axis tight

% led images have no negative values and dont need to be changed in Cscale
if min(data, [], 'all') < 0
    mx = mean(mean(abs(data)));
    mx = double(mx);
    set(gca(), 'CLim', [-1, 1]*mx*5)
end

nMasks = {};
n = true;
while n
    selection = drawfreehand();
    % catch premature ended loop by clicking ESC
    if selection.Position
        msk = createMask(selection);
        nMasks{end+1} = msk;
    else
        n = false;
    end
end
disp(['<> returning cell of ', num2str(n), ' masks'])
