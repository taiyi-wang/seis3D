function [vx_as, dsdt_as, taux_as, tauy_as] = stateODE(Dx_as, s_as, t, p, sig_as, M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% With total displacement up to time 't', compute the shear stress change 
% on the fault. Use force balance based on rate and state friction 
% (aging law) to compute velocity and time derivative of state.

% Input:
% Dx_as = Nz X Nx nodes array of aseismic slips up to current time, t
% s_as  = Nz X Nx nodes array of dimensionless aseismic state variable
% t  = scalar current time
% p  = Nz X Nx nodes array of total pressure on main fault
% sig_as = Nz X Nx normal stress on the main fault
% M   = a structure of relevant constants regarding fault set up (M1)

% Output:
% vx_as   = Nz X Nx array of velocities on the aseismic fault
% dsdt_as = Nz X Nx array of state variable time derivatives (for aseismic main fault)
% taux_as = Nz X Nx array of aseismic shear stress change in x direction 
% tauy_as = Nz X Nx array of aseismic normal stress change in y direction 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute stress change on main fault 
[dtau_as_as_yx, ~, ~] = cmp_traction(Dx_as-t.*M.vpx, [M.mu, M.nu, M.Npad, M.dx, M.dz, M.Nx, M.Nz]); % due to aseismic slip
dtau_as_yx = dtau_as_as_yx; % total shear stress change
tauLockx_as = M.taux_as_0+dtau_as_yx(:,:,1);

% compute velocity in x over whole domain
vx_as = cmp_v(tauLockx_as, s_as, p, sig_as, M);                     % note: increase the normal stress acts against pressure, hence the negative sign
    
% rate of change in state
dsdt_as = (M.b.*M.V0./M.dc).*(exp((M.f0 - s_as)./M.b) - vx_as./M.V0);

% compute total traction in x (shear)
taux_as = tauLockx_as-M.eta.*vx_as; 
    
% compute total traction in y (normal)
tauy_as= sig_as;

end