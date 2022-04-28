function U = eshelby_3dcrack(mode, tau, R, mu, nu, X, Z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the 3D displacement field due to constant stress
% change on a circular crack.

% Input:
% mode   = 1 for antiplane shear in x, 2 for openning mode in y
% tau    = uniform traction change in x
% R      = radius of the circular crack
% mu, nu = shear modulus and Poisson's ratio
% X, Z   = grids of x and z coordinates to compute the displacement at

% Output:
% U      = a structure for displacement field

% Notes:
% For traction in x, a displacement field in x and z is produced.
% For traction in y, a displacement field in y is produced.

% Reference:
% J. D. Eshelby. The determination of the elastic field of an ellipsoidal 
% inclusion, and related problems.Proceedings of the royal society of London. 
% Series A. Mathematical and physical sciences, 241(1226):376–396, 1957.

% To do:
% Check if the z component displacement is correct for the antiplane shear
% case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r2 = X.^2 + Z.^2;
if mode == 1
    tx = tau;
    ux = ((4*(1-nu)*R*tx)/(pi*(2-nu)*mu)).*real(sqrt(1-(r2./(R^2))));
    uz = pi^2 * (1-2*nu) * (2-nu) * R^2 * tx / (16 * mu * (1-nu)).*X;
    U.ux = ux;
    U.uz = uz;
elseif mode == 2
    ty = tau;
    uy = (2*R*(1-nu)*ty/(mu*pi)).*real(sqrt(1 - r2./R^2));
    U.uy = uy;
end

end