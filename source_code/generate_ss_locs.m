function ss_locs = generate_ss_locs(N_ss, sigma_y, xlb, xub, zlb, zub, fault_depth, half_width)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates the spring slider locations randomly distributed in x and z,
% but Gaussian distributed in y (outside of the fault zone)

% Input:
% N_ss = number of spring sliders
% sigma_y = standard deviation (m) for depth (y-dir) distribution of spring sliders
% xlb, xub = lower, upper bound in x (m)
% zlb, zub = lower, upper bound in z (m)
% fault_depth = depth of main fault (positive number, in meters)
% half_width = half width of fault zone (m)

% Output:
% ss_locs = a N_ss X 3 matrix of spring sliders locations (x, z, y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ss_locs = zeros(N_ss, 3);

Lx = xub-xlb;
Lz = zub-zlb;

ss_locs(:,1) = rand(N_ss, 1).*Lx + xlb;
ss_locs(:,2) = rand(N_ss, 1).*Lz + zlb;
dy = normrnd(0,sigma_y, N_ss, 1);
ss_locs(:,3) = sign(dy).* half_width + dy + fault_depth; 

end