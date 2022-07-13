function pPipe = pipeFriction(Q,rho,R,L,f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pressure loss from turbulent pipe friction

% Input:
% Q   = Nt x 1 injection volume history (m^3/s)
% rho = density of fluid 
% R   = well radius
% L   = well length
% f   = Darcy-Weibach frictional coefficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = 2*R;
A = pi*R^2;
v = Q./A;
pPipe = f*L/D*rho.*v.^2./2;
end