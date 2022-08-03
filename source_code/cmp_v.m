
function v = cmp_v(tauLock, s, p, sig_as, M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This funciton computes slip velocity from quasi-dynamic stress balance equation
% Input:
% tauLock = Nz x Nx initial shear stress + loading stress
% s       = Nz x Nx state variable
% p       = Nz x Nx pore pressure
% sig_as  = Nz x Nx normal stress
% M       = structure (M1) containing parameter values for elasticity and friction

% Output:
% v       = Nz x Nx array of velocities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute velocity for aseismic slip
[Nrow, Ncol] = size(tauLock);
v = zeros(Nrow, Ncol);
tauLock_s = abs(tauLock); 

for i = 1:Nrow
    for j = 1:Ncol
        lb = 0; % lower bound of search
        ub = tauLock_s(i,j)/M.eta; % upper bound of search
         
        [V_local, ~] = hybrid(@(vel)solveV(vel,tauLock_s(i, j),s(i,j),M.a(i,j), M.V0, M.eta, sig_as(i,j), p(i,j)),lb,ub, 1e-14, 1e-6);
        
        % just solve scalar velocity
        v(i,j) = V_local.*sign(tauLock(i,j));
    end
end

end