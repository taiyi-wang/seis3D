function Tyx = t_hat1D(kx, kz, Dx, mu, nu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes traction in Fourier domain for 1D slip
% Input:
% kx = wave number in x direction
% kz = wave number in z direction
% Dx = slip in x direction (along fault strike)
% mu = crustal shear modulus
% nu = poisson's ratio

% Output:
% Tyx = a scalar of traction along x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[kx, kz] = meshgrid(kx, kz);
q = sqrt((kx.^2) + (kz.^2));

% Each component Fij represent a component of A*M*A', where A is the
% rotation matrix

F11 = (-mu.*(kx.^2)./(q.*(1-nu))) - ((mu.*(kz.^2))./(q)); % spectral second derivative with regard to x, z
Tyx = F11 .* Dx;

% q = 0 corresponds to kx and kz = 0. In the limit that kx, kz approach zero,
% all the F_ij goes to zero, so corresponding stresses are zero. zero frequency amplitude also corresponds to 
% amplitude of uniform slip over the entire fault surface, which does
% not correspond to stress transfer. 

Tyx(1, 1, :) = 0;
end