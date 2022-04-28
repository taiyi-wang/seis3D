function [x, y] = rad2cat(r, theta, x0, y0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converts coordinates from radial (r, theta) to cartesian (x, y)

% Input:
% r     = radial distance from origin
% theta = azimuth relative to origin
% x0, y0 = x and y coordinates of origin

% OutputL
% x, y   = Nloc X 1 vector of x (or y) locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = r.*cos(theta) + x0;
y = r.*sin(theta) + y0;
end