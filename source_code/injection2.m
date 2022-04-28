function [r, t_s, p_s, pWH, t_s_w0, dpw0, q_s] = injection2(M2, testflag, rglgridflag, plotflag, saveflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform radial linear pressure diffusion simulation given the injection volume
% history at the well (ode15s for solving stiff ODE).

% Input:
% M2       = structure with diffusion problem set up
% testflag    = use step function injection history (1) for testing
% rglgridflag = use regular grid (1) or logrithmic grid (0)
% plotflag = plot in real time simulated well head pressure
% saveflag = save the results in a mat file in the end

% Output:
% r = Npt x 1 radial coordinates
% t_s = Nt x 1 time steps
% p_s = Nt x Npt array of pressure history
% pWh = 1 x Nt well head pressure
% t_s_w0 = Nt x 1(no flow) time steps
% dpw0   = Nt x 1 (no flow) pressure change in well
% q_s = Nt x 1 injection history

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('computing pressure history from injection');

% Load in all parameters---------------------------------------------------
t_max = M2.duration; % simulation time (s)
dt_max = M2.dt_max;     % time step size (s) 

% fault zone properties
w = M2.w;          % width (m)
k = M2.k;          % permeability (m^2)
b = M2.beta_phi;   % porosity*(fluid+pore compressibilities) (1/Pa)
etav = M2.eta_v;   % fluid viscosity (Pa*s)
p0 = M2.p_pore;    % initial pore pressure 

% well properties
Rw = M2.R_w; % radius (m)
Lw = M2.H_w; % length (m)
Sw = 0*pi*Rw^2*Lw*M2.beta_w; % storage capacity (m^3/Pa)
f = M2.f;     % Darcy-Weisbach friction factor
rho = M2.rho; % fluid density (kg/m^3)
g = M2.g;

% unit conversion
Pa2MPa = 1e-6; % convert Pa to MPa
m3s2Ls = 1e3/60; % convert m^3/s to L/min

% Set up spatial discretization (log spacing to capture the near-injector pressure field)---------------
% radial distance, solve for rmin<r<rmax
rmin = Rw;                                                % minimum radius (m)
if testflag == 1
    rmax = M2.rmax;
elseif testflag == 0
    rmax = sqrt(max(abs(M2.X(:))).^2 + max(abs(M2.Z(:))).^2); % maximum radius at least longer than sqrt(xmax^2 + zmax^2)(m) 
end

% spatial discretization (log spacing to capture the near-injector pressure field)
nr = M2.nr; % number of grid intervals = number of grid points - 1
nr = nr + 1;
nr = max(nr,9); % minimum limit for discretization operator

if rglgridflag == 0
    r = (10.^linspace(log10(rmin), log10(rmax), nr))'; % radial coordinate vector (m)
elseif rglgridflag == 1
    h = (rmax-rmin)/(nr-1); % uniform grid spacing (m)
    r = rmin+(0:nr-1)'.*h;
end

% Set up diffusion parameters and boundary conditions-----------------------
B = 2*pi*r.*w.*b;
K = 2*pi*r.*k.*w./etav;
if testflag == 0       % boundary condition at r=rmax
    edgeBC = 'pressure'; 
elseif testflag == 1
    edgeBC = 'flowrate';
end
FZ = Diffusion(nr,B,K,Sw,edgeBC); % diffusion object
FZ = FZ.discretize(r); % spatial discretization

% Set up injection history--------------------------------------------------
if testflag == 1
    nt = ceil(t_max/dt_max);  % number of time steps
    tq = (0:nt)*dt_max;      % time vector (s)
    Q = M2.Q0.*heaviside(tq); 
elseif testflag == 0
    tq = M2.H4_time;
    Q = M2.H4_rate;
end

% Set up stiffness solver-----------------------------------------------------
% calculate the Jacobian for dpdt = f(p1, p2, p3, .. pn, t)
J = FZ.B\(FZ.D2-FZ.Sw.*FZ.S);

dp0 = repelem(0, nr, 1); % initial condition for pressure perturbation
tspan = [0, t_max];      % time span for the solver

options = odeset('RelTol', 1e-2, 'AbsTol', 1e-2, 'MaxStep', dt_max, 'Jacobian', J);
options_noflow = odeset('RelTol', 1e-2, 'AbsTol', 1e-2, 'MaxStep', dt_max);

[t_s, dp] = ode15s(@(t, p) pODE(t, p, FZ, tq, Q), tspan, dp0, options);         % change in pressure field in fault zone
[t_s_w0, dpw0] = ode15s(@(t, p) pODEnoflow(t, p, FZ, tq, Q), tspan, 0, options_noflow); % pressure change in well, if no flow out of well

% wellhead pressure includes pipe friction
q_s = interp1(tq, Q, t_s);
pPipe = pipeFriction(q_s,rho,Rw,Lw,f); % pipe friction pressure loss (Pa)

p_s = dp + p0;                                   % solution at all r in the fault zone
pWH = p_s(:, 1) + pPipe - rho*g*Lw;                    % wellhead pressure with viscous pressure loss(Pa)
    
if plotflag
    figure;
    yyaxis right
    plot(t_s./86400,q_s*m3s2Ls, '-', 'Color', [0.5, 0.5, 0.5]); hold on;
    ylabel('injection rate (m^3/s)');
    yyaxis left
    plot(t_s./86400,pWH.*Pa2MPa, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1); hold on;
    plot(M2.H4_time./86400, M2.H4_wh_p.*Pa2MPa, 'k-', 'LineWidth', 1);
    xlabel('days'); ylabel('well head pressure (MPa)')
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    drawnow
end

if saveflag 
    save('injection_output', 'r', 't_s', 'p_s', 'pWH', 'q_s');     
end

p_s = p_s'; 
t_s = t_s';
q_s = q_s';

disp('finished');
end

function dpdt = pODE(t, p, FZ, tq, Q)
% system of ODE for radial fluid injection problem; spatial discretization
% using SBP
B = FZ.B; S = FZ.S; s = FZ.s; D2 = FZ.D2; Sw = FZ.Sw;

% interpolate injection rate at current time
Qin = interp1(tq, Q, t);

dpdt = B\(s*Qin + (D2-Sw.*S)*p);
end

function dpdt = pODEnoflow(t, p, FZ, tq, Q)
% ODE for the case where the fault zone has zero diffusivity
Sw = FZ.Sw;

% interpolate injection rate at current time
Qin = interp1(tq, Q, t);

dpdt = Qin/Sw;
end