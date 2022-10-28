
function r = cmp_seis_extent(t, dz,  Dtauc, eta, beta,  k, w, Q0, phi, std)
% Semi-analytical solver for spatial extent of fault-weakening resulted 
% seismicity, given the distance 

% Input:
% t  = scalar time 
% dz = scalar vertical distance between the secondary faults and the
%      boundary of main fault
% Dtauc = critical strength drop = f0*effective stress - tau0
% eta   = fluid viscosity
% k, w, phi     = fault zone permeability, width, porosity
% Q0    = scalar value of magnitude of step function injection rate [volume/time]
% std   = standard deviation for pore pressure perturbation away from the
%       fault boundaries

% Output:
% r = maximum (theoretical) radial distance of seismicity caused by direct
% fault weakening

syms r 
Dp = (Q0*eta/(4*pi*k*w))*expint(phi*eta*beta*r^2/(4*k*t))*exp(-dz^2/(2*std^2));

eqn = Dtauc - Dp == 0;
cond = [0, 5e3];

r = double(vpasolve(eqn, cond));

end