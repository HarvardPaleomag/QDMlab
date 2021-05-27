function B = dipole(xc, yc, m, dec, inc, x, y, h)

% rho: radial distance
% phi: polar angle
% z:   height

%% m,de,inc -> mx,my,mz
dec = deg2rad(dec);
inc = deg2rad(inc);

mx = m * cos(dec) * cos(inc);
my = m * sin(dec) * cos(inc);
mz = m * sin(inc);
mvec = [mx, my, mz];

%%
[X, Y, Z] = meshgrid(xc-x, yc-y, h);
% [theta,rho] = cart2pol(X,Y);

B = field(mvec, X, Y, Z);
end

function B = field(m, x, y, z)
r = sqrt(x.^2+y.^2+z.^2);
R = cat(3, x,y,z);
M = ones(size(R));
M(:,:,1) = m(1);
M(:,:,2) = m(2);
M(:,:,3) = m(3);

mu0 = 4*pi*1e-7; %μ0 = 4π×10−7 H/m
a = (mu0 / (4 * pi));
B = a * (3 * dot(M,R,3).*R - r.^2 .* M)./r.^5;
% Bx = (mu0 / (4 * pi)) * (3 * (m(1) * x) .* x - r.^2 * m(1))./r.^5;
% By = (mu0 / (4 * pi)) * (3 * (m(2) * y) .* y - r.^2 * m(2))./r.^5;
% Bz = (mu0 / (4 * pi)) * (3 * (m(3) * z) .* z - r.^2 * m(3))./r.^5;
end
