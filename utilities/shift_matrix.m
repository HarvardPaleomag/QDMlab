function M = shift_matrix(M, dx, dy)
%[M] = shift_matrix(M, dx, dy)
%{
Shifts a mask by dx, dy.
%}
% disp(['<>   shifting Matrix by dx = ' num2str(dx) '; dy = ' num2str(dy)])
if dx <= 0
    dx = abs(dx);
    shifted=zeros(size(M));
    shifted(:,1:end-dx) = M(:, dx+1:end);
elseif dx > 0
    shifted=zeros(size(M));
    shifted(:, dx+1:end) = M(:, 1:end-dx);
end

M = shifted;

if dy <= 0
    dy = abs(dy);
    shifted=zeros(size(M));
    shifted(dy+1:end,:) = M(1:end-dy, :);
elseif dy > 0
    shifted=zeros(size(M));
    shifted(1:end-dy, :) = M(dy+1:end, :);
end

M = shifted;