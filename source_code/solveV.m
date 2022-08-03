
function residual = solveV(V, tauLock,s, a, V0, eta, sigma, p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve for stress-strength residual in the quasi-static stress balance equation
% Input:
% V = scalar slip velocity
% tauLock = scalar initial shear stress + shear stress change due to loading
% a = rate-and-state direct effect coefficient
% V0 = scalar reference velocity
% eta = scalar radiation damping parameter
% sigma = scalar lithostatic normal stress
% p     = scalar pore pressure

% Output:
% residual = scalar stress - strength
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stress = tauLock-eta.*V; % scalar stress magnitude
    
f = a.*asinh(V./(2*V0).*exp(s./a));
strength = f*(sigma-p); 
residual = stress-strength;
end