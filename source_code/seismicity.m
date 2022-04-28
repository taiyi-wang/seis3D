function [t_ss_s, Psi_ss_s, Vx_ss_s, Dx_ss_s, p_ss_s, taux_ss_s, idx_timeout] = seismicity(M1, M2, t_as_s, Dx_as_s, p_as_s, saveflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up randomly located spring sliders driven by pressure diffusion and
% driven by stress transfer due to aseismic slip on the main fault 

% Input:
% M1     = structure containing setup for elasticity and friction
% M2     = structure containing setup for pressure diffusion
% t_as_s = Nt x 1 time steps from adaptive time stepping on main fault
% Dx_as_s = Nz x Nx x Nt array of slip history from main fault
% p_as_s  = Nz x Nx x Nt array of pore pressure history in the fault zone
% save_flag = save the results in a mat file in the end

% Output:
% t_ss_s   = Nss x 1 cell of Ntime x 1 time steps from adaptive time stepping on seismic patches
% Psi_ss_s = Nss x 1 cell of Ntime x 1 dimensionless state variable on each spring slider
% Vx_ss_s  = Nss x 1 cell of Ntime x 1 slip velocity history on each spring slider
% Dx_ss_s  = Nss x 1 cell of Ntime x 1 slip history on each spring slider
% p_ss_s = Nss x 1 cell of Ntime x 1 pore pressure history interpolated to each spring slider
% taux_ss_s = Nss x 1 cell of Ntime x 1 shear stress (along dip) history on each spring slider
% idx_timeout = indices of spring sliders that took too long to compute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('compute off-fault seismicity');

% load 
X = M1.X; Z = M1.Z; 
ss_locs = M1.ss_locs; p_pore = M2.p_pore;
M_as_ss = M1.M_as_ss; sigma_y = M1.sigma_y;
half_width = M1.w/2;
fault_depth = M1.fault_depth;
N_ss = M1.N_ss;
N_as_t = length(t_as_s);

% interpolate the pressure and stress evolution on the main fault to that
% on the spring sliders

dtaux_ss = nan(N_as_t, N_ss); % change in shear stress (along dip) on spring sliders due to main fault slip
dtauy_ss = nan(N_as_t, N_ss); % change in normal stress on spring sliders due to main fault slip
p_ss = nan(N_as_t, N_ss);


for i = 1:N_as_t
    Dx_as = Dx_as_s(:,:,i);
    Dx_as = Dx_as(:);
    dtaux = M_as_ss{1} * Dx_as; 
    dtauy = M_as_ss{2} * Dx_as; 
    
    % interpolate pressure diffusion if not provided, this step can be time
    % consuming
    p = nan(1, N_ss);
    for j = 1:N_ss
        p_fault = interp2(X, Z, squeeze(p_as_s(:, :, i)), ss_locs(j, 1), ss_locs(j, 2));            % interpolate pressure on main fault to the x and z coordinates of the spring sliders
        p(1, j) = interp_pressure(sigma_y, ss_locs(j,3), fault_depth, p_fault, p_pore, half_width); % interpolate pressure on main fault p(x_ss, z_ss, y_as) to p(x_ss, z_ss, y_ss)
    end
    p_ss(i, :) = p';

    dtaux_ss(i, :) = dtaux';
    dtauy_ss(i, :) = dtauy';
end

%% Compute spring slider slip history
% Set up spring slider parameters
% simulation parameters
M.tol = M1.tol;
M.dt = M1.dt;
M.dtmax = M1.dt_max;
M.safety = M1.safetyFactor;

% rate-and-state parameters
M.f0 = M1.f0; % reference friction coefficient
M.V0 = M1.V0; % reference slip velocity
M.a = M1.a_ss; % direct effect parameter
M.b = M1.b_ss; % state evolution parameter
M.dc = M1.dc_ss; % (m) state evolution distance

% elasticity parameters
M.eta = M1.eta; % radiation-damping coefficient (Pa*s/m)
k_x = M1.k_x;     % stiffness of the spring slider

% Initial normal stress
tauy_ss_0 = M1.tauy_ss_0; % (Pa)

% initial conditions
D0 = M1.Dx0_ss; % slip (m)
M.tau0 = M1.taux_ss_0; % shear stress (Pa)
Psi0 = M1.s0_ss; % state variable

% time
tmax = max(t_as_s);

t_ss_s = cell(N_ss, 1);
Dx_ss_s = cell(N_ss, 1);
Psi_ss_s = cell(N_ss, 1);
Vx_ss_s = cell(N_ss, 1);
p_ss_s = cell(N_ss, 1);
taux_ss_s = cell(N_ss, 1);

idx_timeout = []; % record indices of spring sliders that timed out

for i = 1:N_ss
    tStart = tic;
    [ta, Va, Da, Psia, Pa, taua] = springslider(M, D0, Psi0, tauy_ss_0(i), tmax, t_as_s, p_ss(:,i), dtaux_ss(:,i), dtauy_ss(:,i), k_x(i));
    tCmp = toc(tStart); % computational time
    if tCmp >= 40
        idx_timeout = [idx_timeout, i];
    end
    t_ss_s{i} = ta;
    Dx_ss_s{i} = Da;
    Psi_ss_s{i} = Psia;
    Vx_ss_s{i} = Va;
    p_ss_s{i}= Pa;
    taux_ss_s{i} = taua;
end

%% Store results
if saveflag 
    save('ss_output', 't_ss_s', 'Dx_ss_s', 'Psi_ss_s', 'Vx_ss_s', 'p_ss_s', 'taux_ss_s', 'idx_timeout', 'dtaux_ss', 'dtauy_ss', 'p_ss');
end

disp('finished')
end














