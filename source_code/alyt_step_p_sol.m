function p_t = alyt_step_p_sol(time, r0, Q0, k, beta, phi, eta_v, w, plot01)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytical solution to step function injection problem.
% Solution in polar coordinatesis obtained by convolving the Green's function to the 2D
% diffusion problem with a step function injection history

% Reference:
% Earlougher RC Jr. Advances in well test analysis. SPE Monograph. Vol.5,
% SPE of AIME, Dallas, 1977

% Taiyi Wang 07/14/21

% Input:
% time = a Ntime X 1 vector of time steps
% r0   = radius at which to compute pressure history
% Q0   = injection rate (m^3/s)
% k    = fault zone permeability (m^2)
% beta = compressibility (Pa^-1)
% phi  = porosity
% eta_v = fluid viscosity (Pa s)
% w     = damage zone width
% plot01 = 1 for plotting 

% Output:
% p_t = predicteed pressure history
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = @(r,t) Q0.*eta_v./(4.*pi.*k.*w).*(-ei(-(phi.*eta_v.*beta.*r.^2)./(k.*t)./4));
p_t = p(r0, time);


if plot01 == 1
    plot(time, p_t)
end

end