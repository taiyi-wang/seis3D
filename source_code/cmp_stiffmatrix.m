function [M_ss_as, M_as_ss, M_ss_ss] = cmp_stiffmatrix(as_locs, ss_locs, as_dx, as_dz, ss_dx, ss_dz, mu ,nu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the stiffness matrices relating slip on main fault or
% spring-sliders with traction changes on each of the faults

% Input:
% as_locs = N X 3 array of coordinates [x, y, z] for aseismic slip grid points (all input depth should be positive)
% ss_locs = M X 3 array of coordinates [x, y, z] for seismic slip grid points (all input depth should be positive)
% as_dx, as_dz = x, z direction dimension of main fault elements
% ss_dx, ss_dz = x, z direction dimension of spring slider sizes
% mu = shear modulus
% nu = Poisson's ratio

% Output:
% M_ss_as = an array of 2 M X N matrices relating slip on spring sliders to shear and normal tractions changes on the main fault
% M_as_ss = an array of 2 N X M matrices relating slip on the main fault to shear and normal tractions changes on the spring sliders
% M_ss_ss = an array of 2 M X M matrices relating slip on the spring sliders to shear and normal tractions changes on the spring sliders

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N_as, ~] = size(as_locs);
[N_ss, ~] = size(ss_locs);

% 1. Form the matrix linking seismic slip to aseismic stress change
M_ss_as_yx = zeros(N_as, N_ss);
M_ss_as_yy = zeros(N_as, N_ss);
%{
for i = 1:N_ss
    ss_loc_x = ss_locs(i, 1)+ss_dx/2; % source location is defined at the right boundary of the displocation
    ss_loc_z = ss_locs(i, 2);
    ss_loc_y = ss_locs(i, 3); % positive source depth is "in the ground"
    
    as_loc_x = as_locs(:, 1);
    as_loc_z = as_locs(:, 2);
    as_loc_y = -as_locs(:, 3); % negative receiver depth is "in the ground"
    
    strike_slip = 0;
    dip_slip = -1; % unit slip [m] towards positive-x
    openning = 0; 
    
    dip = 0; strike = 0; 
    mdl_ss = [ss_dx, ss_dz, ss_loc_y, dip, strike, ss_loc_x, ss_loc_z, strike_slip, dip_slip, openning]'; % source parameters
    [~, ~, S, ~] = disloc3dpm(mdl_ss, [as_loc_x, as_loc_z, as_loc_y]', mu, nu);
    
    % Extract stress 
    M_ss_as_yx(:, i) = transpose(S(3,:)); % sigma_yx
    M_ss_as_yy(:, i) = transpose(S(6,:)); % sigma_yy
end
%}
M_ss_as = {M_ss_as_yx, M_ss_as_yy};

disp('finish setting up stiffness matrix M_ss_as')

% 2. Form the matrix linking aseismic slip to seismic stress change
M_as_ss_yx = zeros(N_ss, N_as);
M_as_ss_yy = zeros(N_ss, N_as);

parfor j = 1:N_as
    fprintf('at %f percent;', j/N_as*100);
    
    as_loc_x = as_locs(j, 1)+ as_dx/2; % source location is defined at the right boundary of the displocation
    as_loc_z = as_locs(j, 2);
    as_loc_y = as_locs(j, 3); % positive source depth is "in the ground"
    
    ss_loc_x = ss_locs(:, 1);
    ss_loc_z = ss_locs(:, 2);
    ss_loc_y = -ss_locs(:, 3); % negative receiver depth is "in the ground"
    
    strike_slip = 0;
    dip_slip = -1; % unit slip [m] towards positive-x
    openning = 0; 
    
    dip = 0; strike = 0; 
    mdl_as = [as_dx, as_dz, as_loc_y, dip, strike, as_loc_x, as_loc_z, strike_slip, dip_slip, openning]';
    [~, ~, S, ~] = disloc3dpm(mdl_as, [ss_loc_x, ss_loc_z, ss_loc_y]', mu, nu);
    
    % Extract stress 
    M_as_ss_yx(:, j) = transpose(S(3,:));
    M_as_ss_yy(:, j) = transpose(S(6,:));
end

M_as_ss = {M_as_ss_yx, M_as_ss_yy};

disp('finish setting up stiffness matrix M_as_ss')
disp('start setting up stiffness matrix M_ss_ss')

% 3. Form the matrix linking seismic slip to stress change at seismic patch
M_ss_ss_yx = zeros(N_ss, N_ss);
M_ss_ss_yy = zeros(N_ss, N_ss);

parfor k = 1:N_ss
    % seismic slip source
    ss_loc_x_s = ss_locs(k, 1) + ss_dx/2; % source location is defined at the right boundary of the displocation
    ss_loc_z_s = ss_locs(k, 2);
    ss_loc_y_s = ss_locs(k, 3); % positive source depth is "in the ground"
    
    % seismic slip receiver (put at the center of the seismic fault patch)
    ss_loc_x_r = ss_locs(:, 1); 
    ss_loc_z_r = ss_locs(:, 2);
    ss_loc_y_r = -ss_locs(:, 3); % negative receiver depth is "in the ground"
    
    strike_slip = 0;
    dip_slip = -1; % unit slip [m] towards positive-x
    openning = 0; 
    
    dip = 0; strike = 0; 
    mdl_as = [ss_dx, ss_dz, ss_loc_y_s, dip, strike, ss_loc_x_s, ss_loc_z_s, strike_slip, dip_slip, openning]';
    [~, ~, S, ~] = disloc3dpm(mdl_as, [ss_loc_x_r, ss_loc_z_r, ss_loc_y_r]', mu, nu);
    
    % Extract stress 
    M_ss_ss_yx(:, k) = transpose(S(3,:));
    M_ss_ss_yy(:, k) = transpose(S(6,:));
end

M_ss_ss = {M_ss_ss_yx, M_ss_ss_yy};

disp('finish setting up stiffness matrix M_ss_ss')

end
