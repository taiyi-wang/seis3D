function field_interp = interp_t(t, field, tq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolate scalar field in time
% Input:
% t     = Nt x 1 vector of time corresponding to 'field'
% field = Nx x Nz x Nt array
% tq    = scalar query time at which the scalar field is sampled

% Output:
% field_interp = Nx x Nz array interpolated at 't'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Nz, Nx, ~] = size(field);
field_interp = nan(Nz, Nx);

for i = 1:Nz
    for j = 1:Nx
        field_interp(i, j) = interp1(t, squeeze(field(i, j, :)), tq);
    end
end

end