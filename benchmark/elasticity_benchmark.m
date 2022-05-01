% Check that elasticity is working properly with 1. spectral element method
% and 2. Okada discrete elements
%% Add paths 
current_dir  = pwd;
idcs   = strfind(current_dir ,'/');
above_dir = current_dir(1:idcs(end)-1); % one level above current director

addpath(fullfile(above_dir, '/source_code'));

%% 1. Compare with analytical Eshelby solution
% Compute Eshelby constant stress drop analytical displacement field, feed
% that into the 3D stress solver, cmp_traction.m, and check if a constant
% stress drop field is reproduced.

% Note here I use the parameter and domain size consistent with Cooper
% Basin simulation to demonstrate that the mesh resolution used is
% approproate

mode = 1; % 1 for antiplane shear, 2 for openning mode

mu = 2.4e10; nu = 0.25; Npad = 0;

% Dimension of model domain (m)
xL = -10e3; xU = 10e3; 
zL = -10e3; zU = 10e3; 

Lx = xU - xL;
Lz = zU - zL;

R = 5e3; % dimension of the crack radius 

% Number of nodes on the physical grid
Nx = 215; Nz = 215;

% Define fault patch
dx = Lx/(Nx-1);
dz = Lz/(Nz-1);
x_locs = linspace(xL,xU, Nx);
z_locs = linspace(zL,zU, Nz);
[X, Z] = meshgrid(z_locs, x_locs);

consts = [mu, nu, Npad, dx, dz, Nx, Nz];

% Compute displacement distribution from constant stress drop circular
% crack
U1 = eshelby_3dcrack(mode, 1e6, R, mu, nu, X, Z);
ux_1 = U1.ux;

[T_field, kx, kz] = cmp_traction(ux_1, consts);

figure;
subplot(1, 2, 1)
surf(X./1e3, Z./1e3, ux_1, 'edgecolor', 'none');colorbar;
view([0, 90]); xlabel('x (km)'); ylabel('z (km)');
title('\delta_x (m)'); xlim([xL./1e3, xU./1e3]);ylim([zL./1e3, zU./1e3]);

subplot(1, 2, 2);
surf(X./1e3, Z./1e3, T_field(:,:,1)./1e6, 'edgecolor', 'none');colorbar
view([0, 90]); xlabel('x (km)'); ylabel('z (km)');
title('\Delta \tau_{yx} (MPa)'); xlim([xL./1e3, xU./1e3]);ylim([zL./1e3, zU./1e3]);

%suptitle(fprintf('%i by %i nodes used for numerical calculation', Nz, Nx))

%% Compute error as a function of number of nodes in each direction and zero padding
% demonstrate convergence with mesh refinement and zero padding
N_list = 11:16:301;                         % have to be odd numbers
pad_L = 10e3;                                 % zero-padd this much distance on either side of domain
err_array = cell(length(N_list), 1);   
Npad_array = cell(length(N_list), 1);
cross_sections = cell(length(N_list), 1);      
alpha_list = linspace(0.1, 1, length(N_list)); % a list of different levels of transparency for plotting

for i = 1:length(N_list)
    Nz = N_list(i);
    Nx = N_list(i);
    dx = Lx/(Nx-1);
    dz = Lz/(Nz-1);
    x_locs = linspace(xL,xU, Nx);
    z_locs = linspace(zL,zU, Nz);
    [X, Z] = meshgrid(z_locs, x_locs);
    U1 = eshelby_3dcrack(mode, 1e6, R, mu, nu, X, Z); % analytical result
    ux_1 = U1.ux;
    Npad_list = 0:1:ceil(pad_L/dx);
    Npad_array{i} = Npad_list;
    err_for_Npad = zeros(1, length(Npad_list)); % temporarily store errors for this mesh resolution, at various numbers of zero padding
    for j = 1:length(Npad_list)
        Npad = Npad_list(j);
        consts = [mu, nu, Npad, dx, dz, Nx, Nz];
        [T_field, ~, ~] = cmp_traction(ux_1, consts); % Numerical result
        T_inside = T_field(sqrt(X.^2+Z.^2)<R-100); % traction calculations within the circular crack
        err_for_Npad(j) =  norm(T_inside - (-1e6))/sqrt(Nz*Nx);
    end
    err_array{i, 1} = err_for_Npad;
    cross_sections{i} = T_field((Nz+1)/2, :);
end

figure;
syms x
fplot(-rectangularPulse(-R, R, x), 'r', 'LineWidth', 1); hold on;
for i = 1:length(N_list)
    plot1 = plot(linspace(xL./1000, xU./1000, N_list(i)), cross_sections{i}./1e6, 'k', 'LineWidth', 1);
    plot1.Color(4) = alpha_list(i);
end
ylabel('shear stress change (MPa)')
xlim([xL./1000, xU./1000])
xlabel('x (km)');


figure;
for i = 1:length(N_list)
    Npad_list = Npad_array{i};
    ratios = (2*Npad_list.*Lx/(N_list(i)-1)+Lx)./Lx; % ratio of zero-padded domain length over real domain length
    plot2 = plot(ratios, err_array{i}./1e6, 'k-', 'LineWidth', 1); 
    plot2.Color(4) = alpha_list(i);hold on;
    xlabel('ratio of zero-padded domain length over real domain length');ylabel('RMS error in stress change (MPa)');
end

%% 2. Compute the displacement and traction field induced by a rectangular Okada element

% Show that the right sense of slip is used at the source of the Okada code
% Here I compute the displacements and stress change at two depths
% locations above and below the source, respectively.

% Note observation depth should be set to be negative (within the earth),
% whereas the source depth should be set to positive to indicate within the
% earth
mu = 2.4e10;   
nu = 0.25;

% Define spring slider
ss_dx = 20; ss_dz = 20; ss_dy = 20;  % dimension
ss_loc_x = 0; ss_loc_z = 0; ss_loc_y = 4050; % location
dip = 0; strike = 0; dip_slip = -1; openning =0; strike_slip = 0; % sense of motion

mdl_ss = [ss_dx, ss_dz, ss_loc_y, dip, strike, ss_loc_x+ss_dx/2, ss_loc_z, strike_slip, dip_slip, openning]'; % source parameters

% Define receiver locations

xL = -1e2; xU = 1e2; 
zL = -1e2; zU = 1e2; 
yL = ss_loc_y+100; yU = ss_loc_y-100; % depth 

Lx = abs(xU - xL);
Lz = abs(zU - zL);
Ly = abs(yU - yL);

Nx = 70;
Nz = 70;
Ny = 70;

re_dx = Lx/(Nx-1);
re_dz = Lz/(Nz-1);

x_locs = linspace(xL,xU, Nx);
y_locs = linspace(yL,yU, Ny);
z_locs = linspace(zL,zU, Nz);

[X, Z] = meshgrid(x_locs, z_locs);
[~, Y] = meshgrid(x_locs, y_locs);
   
re_locs_l = [X(:), Z(:), repelem(-yL, Nx*Nz)']; % below spring slider
re_locs_u = [X(:), Z(:), repelem(-yU, Nx*Nz)']; % above spring slider
re_locs_v = [X(:), repelem(0, Nx*Ny)', -Y(:)];     % vertical plane (x-y) across spring slider

% compute displacements at x-z planes below and above spring slider
[U_l, ~, S_l, ~] = disloc3dpm(mdl_ss, re_locs_l', mu, nu);
[U_u, ~, S_u, ~] = disloc3dpm(mdl_ss, re_locs_u', mu, nu);

% compute displacements at the x-y plane
[U_v, ~, S_v, ~] = disloc3dpm(mdl_ss, re_locs_v', mu, nu);
      
figure;
quiver3(X, Z, -yL.*ones(Nz, Nx), reshape(U_l(1,:), [Nz, Nx]), reshape(U_l(2,:), [Nz, Nx]), reshape(U_l(3,:), [Nz, Nx]), 'b', 'LineWidth', 1); hold on;
quiver3(X, Z, -yU.*ones(Nz, Nx), reshape(U_u(1,:), [Nz, Nx]), reshape(U_u(2,:), [Nz, Nx]), reshape(U_u(3,:), [Nz, Nx]), 'r', 'LineWidth', 1);
xlabel('x (m)'); ylabel('z (m)'); zlabel('y (m)'); view([0, 90]);
title('displacement vectors on planes 100 meters away from the main fault')

% compute traction changes on planes above, below, and vertically transect
% the source
re_locs_l(:,3) = abs(re_locs_l(:,3));
[~, M_re_l, ~] = cmp_stiffmatrix(re_locs_l, [ss_loc_x, ss_loc_z, ss_loc_y], re_dx, re_dz, ss_dx, ss_dz, mu ,nu);

re_locs_u(:,3) = abs(re_locs_u(:,3));
[~, M_re_u, ~] = cmp_stiffmatrix(re_locs_u, [ss_loc_x, ss_loc_z, ss_loc_y], re_dx, re_dz, ss_dx, ss_dz, mu ,nu);

re_locs_v(:,3) = abs(re_locs_v(:,3));
[~, M_re_v, ~] = cmp_stiffmatrix(re_locs_v, [ss_loc_x, ss_loc_z, ss_loc_y], re_dx, re_dz, ss_dx, ss_dz, mu ,nu);

figure; 

subplot(2,3,1)
surf(X, Z, reshape(M_re_u{1}./1e6, Nz, Nx), 'EdgeColor', 'none'); view([0, 90]); colorbar;
hold on; plot(0, 0, 'k*');xlabel('x (m)'); ylabel('z (m)');
title('\Delta \tau_{yx} (MPa) above source');

subplot(2,3,4)
surf(X, Z, reshape(M_re_l{1}./1e6, Nz, Nx), 'EdgeColor', 'none'); view([0, 90]); colorbar;
hold on; plot(0, 0, 'k*');xlabel('x (m)'); ylabel('z (m)');
title('\Delta \tau_{yx} (MPa) on x-z plane below source');

subplot(2,3,2)
surf(X, Z, reshape(M_re_u{2}./1e6, Nz, Nx), 'EdgeColor', 'none'); view([0, 90]); colorbar;
hold on; plot(0, 0, 'k*');xlabel('x (m)'); ylabel('z (m)');
title('\Delta \tau_{yy} (MPa) above source');

subplot(2,3,5)
surf(X, Z, reshape(M_re_l{2}./1e6, Nz, Nx), 'EdgeColor', 'none'); view([0, 90]); colorbar;
hold on; plot(0, 0, 'k*');xlabel('x (m)'); ylabel('z (m)');
title('\Delta \tau_{yy} (MPa) on x-z plane below source');

subplot(2,3,3)
surf(X, -Y, reshape(M_re_v{1}./1e6, Nz, Nx), 'EdgeColor', 'none'); view([0, 90]); colorbar;
hold on; plot(0, 0, 'k*');xlabel('x (m)'); ylabel('y (m)');ylim([-yL, -yU])
title('\Delta \tau_{yx} (MPa) on x-y plane');

subplot(2,3,6)
surf(X, -Y, reshape(M_re_v{2}./1e6, Nz, Nx), 'EdgeColor', 'none'); view([0, 90]); colorbar;
hold on; plot(0, 0, 'k*');xlabel('x (m)'); ylabel('y (m)');ylim([-yL, -yU])
title('\Delta \tau_{yy} (MPa) on x-y plane');



