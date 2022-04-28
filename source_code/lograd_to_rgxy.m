function p_s_xy = lograd_to_rgxy(r, t_s, p_s, M1, saveflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert stored axially-symmetric field onto a regular, rectangular grid

% Input:
% r   = Npt x 1 radial coordinates
% t_s = 1 x Nt time steps
% p_s = Npt x Nt array of pressure history
% M1  = structure containing set up for elasticity and friction
% saveflag = 1 for saving output

% Output:
% p_s_xy = Nz x Nx x Nt array for pressure evolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('converting pressure history from radial to Cartesian coordinates');

X = M1.X; Z = M1.Z;
R = sqrt(X.^2 + Z.^2);  % radial distance from injector to grid points
Nt = length(t_s);       % number of time steps
Nz = M1.Nz; Nx = M1.Nx;

p_s_xy = nan(Nz, Nx, Nt);
for i = 1:Nt
    pr = p_s(:, i);
    for m = 1:Nz
        for n = 1:Nx
            if R(m, n) > 0
                p_s_xy(m, n, i) = interp1(r, pr, R(m, n));
            elseif R(m, n) == 0
                % assume pressure at R = 0 is the same as that at R = well
                % radius 
                p_s_xy(m, n, i) = pr(1);
            end
        end
    end
end

if saveflag == 1
    save('p_s_xy.mat', 'p_s_xy')
end

disp('finished');
end