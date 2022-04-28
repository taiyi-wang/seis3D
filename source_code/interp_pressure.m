function p = interp_pressure(sigma_y, ss_depth, fault_depth, p_fault, p0, half_width)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function interpolates pressure off fault according to a given normal
% distribution

% Input:
% sigma_y     = standard deviation (m) for pressure (y-dir) distribution 
% ss_depth    = depths (m) of spring slider(s)
% fault_depth = depth (m) of main fault
% p_fault     = pressure(s) (Pa) at main fault location closest to spring slider
% p0          = initial uniform pressure on the main fault (Pa)
% half_width  = half width of fault zone

% Output:
% p = pressure inferred at the depth of spring slider (Pa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dp = p_fault - p0;

dy = abs(ss_depth-fault_depth);
p = dp.*exp(-0.5*((dy - half_width)/sigma_y).^2) + p0;
end