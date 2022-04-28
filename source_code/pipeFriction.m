function [pPipe,tau] = pipeFriction(Q,rho,R,L,f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pressure loss from turbulent pipe friction

% Input:
% Q   = Nt x 1 injection volume history (m^3/s)
% rho = density of fluid 
% R   = well radius
% L   = well length
% f   = Darcy-Weibach frictional coefficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = pi*R^2; % cross-sectional area
tau = f/4*rho*(Q/A).^2; % wall shear stress
pPipe = tau*2*L/R; % pressure loss (pi*R^2*p/L=2*pi*R*tau)
end