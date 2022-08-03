function [M1, M2] = setup_model()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up parameters for simulations with default parameters

% Output:
% M1 = a structure for rate-and-state fault parameters
% M2 = a structure for pressure diffusion and well parameters

% Taiyi Wang 09/03/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% General parameters
duration = 17*86400; % 17 days of injection

Npad = 0;            % number of grid points added for each side of zero-padding for computing traction
dt_max = 86400/24/2; % maximum time step size
dt_min = 0.1;        % minumum time step size ~ radiation damping parameter/fault stiffness
safetyFactor = 0.8;  % safty factor for adaptive time stepping
tol = 10^(-6);       % tolerance for time stepping
dt = 1;              % initial time step size
%% 1. Set up elasticity/friction parameters 

% Rate and state
f0 = 0.6;  % reference friction coefficient
V0 = 1e-6; % reference slip velocity
dc = 1.53e-5; % state evolution distance [m] (Marone, 1998)

% Elasticity
mu = 2.4e10;   % shear modulus of granite [Pa] (Hoek and Bray, 1981)
nu = 0.25;     % Poisson's ratio
eta = 4550000; % Radiation damping parameter

% Dimension of model domain (reference is Fig. 2a of Baisch 2019) 
xL = -10e3; xU = 10e3; % along dip, 10 km default(6 degree towards west)
zL = -10e3; zU = 10e3; % along strike, 10 km default (204 degrees)

Lx = xU - xL;
Lz = zU - zL;

% Number of nodes on the physical grid (main fault)
Nx = 215; 
Nz = 215;
as_y = 4100; % depth of main fault = abs(depth)

% Define fault patch
dx = Lx/(Nx-1);
dz = Lz/(Nz-1);
x_locs = linspace(xL,xU, Nx);
z_locs = linspace(zL,zU, Nz);
[X, Z] = meshgrid(x_locs, z_locs);

% Spatially variable parameters
a = repelem(0.015, Nz, Nx);
b = repelem(0.012, Nz, Nx);

% tectonic loading velocity                               
vpx = 0;                                 

% Initial state variable
s0 = repelem(0.75, Nz, Nx);    

% Initial slip
Dx0 = zeros(Nz, Nx);                               % total slip in x

%% 2. Setup stiffness matrix for one-way coupled aseismic->seismic slip
% spatial distribution of spring-slider locations
N_ss = 1000;  % number of spring sliders
sigma_y = 2; % standard deviation in y for spring slider distribution (m)
xL = -1e3; xU = 1e3; % along dip (6 degree towards west)
zL = -1e3; zU = 1e3; % along strike (204 degrees)
w = 6;         % width of fault zone (m) (Holl and Barton, 2015)

ss_locs = generate_ss_locs(N_ss, sigma_y, xL, xU, zL, zU, as_y, w/2);

% frictional properties
a_ss = 0.015; % direct effect parameter
b_ss = 0.018; % state evolution parameter
dc_ss = 1.53e-5; % state evolution distance (m) 
s0_ss = 0.75;  % initial state variable

% Compute the stiffness matrices relating slip on main fault to traction
% changes at spring sliders
as_locs = [X(:), Z(:), repelem(as_y, Nx*Nz)']; % location of the aseismic patches
as_dx = dx; as_dz = dz;       % size of aseismic slip patch
ss_dx = 20; ss_dz = 20;       % size of seismic slip patch

%[~, M_as_ss, M_ss_ss] = cmp_stiffmatrix(as_locs, ss_locs, as_dx, as_dz, ss_dx, ss_dz, mu, nu);
%% 3. Pressure diffusion/permeability evolution parameters
nr = 700; % number of grid intervals for the radial pressure diffusion solver (at least 500 for R = 10 km to ensure convergence)

% Habanero 4 (2012)
R_w = 0.1778/2;     % radius of tubing [m] (well completion report, section 4.0, Pg. 23)
H_w = 4077;         % well depth [m] (Holl & Barton 2015, Fig. 2)
V_w = pi*R_w^2*H_w; % volume of well [m^3]
k = 7e-13;          % permeability [m^2] (from trial and error)
f = 0.01;           % Darcy-Weishbach friction factor (from trial and error)
beta = 1e-8;        % Fluid plus pore compressibility [Pa^-1]
phi = 1e-2;         % Pore volume fraction
beta_w = 4.41e-10;  % compressibility of well [Pa] (assume water dominates compressibility, at 10 MPa, and 25 degree Celsius, Fine & Millero, 1973)
eta_v = 8.90e-4;    % Fluid viscosity [Pa s] (assume water, at 25 degree celsius)
beta_phi = beta*phi; % compressibility X pore volume fraction (the two parameters only appear as product in the governing equations)
rho = 1e3;              % fluid density [kg m^3] (assume water)
g = 9.8;                % gravitaional acceleration m/s^2

% coordinates of the injection well
x0_4 = 0; z0_4 = 0;       % Habanero 4

% Resolve the stress on the fault, which is at a median depth of ~4100m
% (values from Holl and Barton, 2015; Fig. 6)
theta = 10/180*pi;              % fault dip is 10 degrees
s_v = 100e6;                    % lithostatic stress [Pa]
p_pore = 73.82e6;               % hydrostatic pressure + over pressure [Pa]
s_Hmax = 160e6;                 % maximum horizontal stress [Pa] 

rho_r = 2700;                % bulk density of granite [kg m^-3]

[~, snt,snn]=rotate_stress(-s_Hmax,0,-s_v,theta);
taux_as_0 = 15e6;                        % Initial shear stress on the main fault [Pa]
tauy_as_0 = repelem(abs(snn), Nz, Nx);   % initial uniform main fault normal stress [Pa]
taux_ss_0 = 15e6;                         % Initial shear stress on the spring slider [Pa]
tauy_ss_0 = s_v-rho_r*g*(as_y - ss_locs(:, 3));               % initial spring slider normal stress [Pa]
%%  Set parameters to user specified values
% M1 is the structure for parameters relevant for the elasticity and friction problem

% Elasticity/friction
M1.duration = duration;
M1.Npad = Npad;
M1.safetyFactor = safetyFactor;
M1.tol = tol;
M1.dt = dt;
M1.dt_max = dt_max;
M1.dt_min = dt_min;
M1.N_ss = N_ss;

M1.f0 = f0;
M1.V0 = V0; 
M1.dc = dc; 

M1.mu = mu; 
M1.nu = nu; 
M1.eta = eta; 
M1.tauy_as_0 = reshape(tauy_as_0, [Nz, Nx]); 
M1.tauy_ss_0 = repelem(tauy_ss_0, N_ss)';

M1.a = a; M1.b = b;
M1.Lx = Lx; M1.Lz = Lz; M1.dx = dx; M1.dz = dz;M1.X = X; M1.Z = Z; M1.fault_depth = as_y;
M1.Nx = Nx; M1.Nz = Nz; 
M1.taux_as_0 = taux_as_0; 
M1.tauy_as_0 = tauy_as_0; 
M1.s0 = s0;
M1.Dx0 = Dx0;
M1.k_x = diag(M_ss_ss{1});

M1.vpx = vpx; 

M1.w = w; M1.sigma_y = sigma_y;

% spring sliders
M1.ss_locs = ss_locs;
M1.a_ss = a_ss; M1.b_ss = b_ss;
M1.dc_ss = dc_ss;
M1.M_as_ss = M_as_ss;
M1.taux_ss_0 = taux_ss_0;
M1.tauy_ss_0 = tauy_ss_0;
M1.ss_dx = ss_dx; M1.ss_dz = ss_dz;
M1.Dx0_ss = 0; % initial slip 
M1.s0_ss = s0_ss; % initial state variable 

% M2 is the structure for all parameters relevant for the diffusion problem
M2.duration = duration;
M2.dt_max = dt_max;
M2.nr = nr;

M2.X = X; M2.Z = Z; M2.dx = dx; M2.dz = dz; 

M2.x0_4 = x0_4; M2.z0_4 = z0_4;

% initial conditions
M2.tauy_as_0 = tauy_as_0; 
M2.p_pore = p_pore;

% constants
M2.R_w = R_w; M2.H_w = H_w; M2.V_w =V_w; M2.beta_w = beta_w; M2.f = f;
M2.eta_v = eta_v; M2.beta = beta; M2.phi = phi; M2.beta_phi = beta_phi; M2.k = k;      
M2.rho = rho; M2.w = w; M2.g = g; 

end