function [T_field, kx, kz] = cmp_traction(Sx, consts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine computes the 3D traction on a fault using spectral method

% Input:
% Sx     = 2D array of slip in x directions
% consts = [mu, nu, Npad(number of zero padding on either side of fault), dx(dimension of each fault patch), dz, Nx (number of nodes in each direction), Nz]

% Output
% T_field = traction field produced by the prescribed displacement field
% kx, kz  = wave numbers in x and z directions

% Taiyi Wang 03/09/2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = consts(1);
nu = consts(2);
Npad = consts(3);
dx = consts(4); dz = consts(5);
Nx = consts(6); Nz = consts(7);
Nfx = Nx-1; Nfz = Nz-1;

% zero-pad X, Z, Sx to increase frequency resolution
Nfx_pad = Nfx+2*Npad; Nfz_pad = Nfz+2*Npad;

Sx_pad = zeros(Nfz_pad, Nfx_pad);
Sx_pad(Npad+1:Npad+Nfz, Npad+1:Npad+Nfx) = Sx(1:Nfz, 1:Nfx);

kx_pad = [0:Nfx_pad/2, (-Nfx_pad/2+1):-1]/Nfx_pad*2*pi/dx; % wavenumber
kz_pad = [0:Nfz_pad/2, (-Nfz_pad/2+1):-1]/Nfz_pad*2*pi/dz;

% Compute slip in Fourier domain
Sx_h_pad = fft2(Sx_pad);
T_field_pad = t_hat1D(kx_pad, kz_pad, Sx_h_pad, mu, nu);

% Get stress in spatial domain
T_field_pad = ifft2(T_field_pad, 'symmetric');

% Remove zero padding
T_field = T_field_pad(Npad+1:Npad+Nfz, Npad+1:Npad+Nfx);
kx = kx_pad(Npad+1:Npad+Nfx);
kz = kz_pad(Npad+1:Npad+Nfz);

T_field(end+1, :) = T_field(1,:);
T_field(:, end+1) = T_field(:,1);
end



