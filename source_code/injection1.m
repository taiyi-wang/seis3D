function [r, t_s, p_s, pWH, dpw0, q_s] = injection1(M2, testflag, rglgridflag, plotflag, saveflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform radial linear pressure diffusion simulation given the injection volume
% history at the well (backward Euler for solving ODE).

% Input:
% M2          = structure with diffusion problem set up
% testflag    = use step function injection history (1) for testing
% rglgridflag = use regular grid (1) or logrithmic grid (0)
% plotflag = plot in real time simulated well head pressure
% saveflag = save the results in a mat file in the end

% Output:
% r = Npt x 1 radial coordinates
% t_s = 1 x Nt time steps
% p_s = Npt x Nt array of pressure history
% pWh = 1 x Nt well head pressure
% q_s = 1 x Nt injection history

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in all parameters
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

% radial coordinate vector (m)
if rglgridflag == 0
    r = (10.^linspace(log10(rmin), log10(rmax), nr))'; % logirithmic grid spacing (m)
elseif rglgridflag == 1
    h = (rmax-rmin)/(nr-1); % uniform grid spacing (m)
    r = rmin+(0:nr-1)'.*h;
end

% spatial discretization and diffusion parameters
B = 2*pi*r.*w.*b;
K = 2*pi*r.*k.*w./etav;
if testflag == 0       % boundary condition at r=rmax
    edgeBC = 'pressure'; 
elseif testflag == 1
    edgeBC = 'flowrate';
end
FZ = Diffusion(nr,B,K,Sw,edgeBC); % diffusion object
FZ = FZ.discretize(r); % spatial discretization

dp = repelem(0, nr, 1); % initial condition for pressure perturbation

% time stepping
tmax = M2.duration; % simulation time (s)
dt = M2.dt_max;     % time step size (s)
nt = ceil(tmax/dt);   % number of time steps
t_s = (0:nt)*dt;      % time vector (s)

% injection history volume rate (m^3/s)
if testflag == 1
    Q = M2.Q0.*heaviside(t_s); 
elseif testflag == 0
    Q = interp1(M2.H4_time, M2.H4_rate, t_s);
end
q_s = Q; 

% time step, optionally saving and plotting solution
Pa2MPa = 1e-6; % convert Pa to MPa
m3s2Ls = 1e3/60; % convert m^3/s to L/min

% arrays storing solution at all time steps
p_s = nan(nr,nt+1); p_s(:,1) = dp + p0;                       % solution at all r
dpw0 = nan(1,nt+1); dpw0(1) = 0; % pressure change in well, if no flow out of well
% wellhead pressure includes pipe friction
pPipe = pipeFriction(Q,rho,Rw,Lw,f); % pipe friction pressure loss (Pa)

if plotflag, figure(1),clf, end

for n=1:nt
    % backward Euler update, note that Q(n) is at t+dt
    dp = FZ.update_diffusion(dp,dt,Q(n));
    p_s(:,n+1) = dp + p0;
    pWH = p_s(1, 1:n+1) + pPipe(1:(n+1))- rho*g*Lw; % wellhead pressure (Pa)
    dpw0(n+1) = dpw0(n)+dt*Q(n+1)/Sw;               % pressure change in the well without outflow
    
    if plotflag
        yyaxis right
        plot(t_s./86400,Q*m3s2Ls, '-', 'Color', [0.5, 0.5, 0.5]); hold on;
        ylabel('injection rate (L/s)');
        yyaxis left
        plot(t_s(1:(n+1))./86400,pWH*Pa2MPa, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1); hold on;
        %plot(t_s(1:(n+1))./86400,p_s(1, 1:n+1)*Pa2MPa, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1); hold on;
        %plot(t_s./86400,dpw0*Pa2MPa, '-', 'Color', [0.9856, 0.7372, 0.2537], 'LineWidth', 1)
        plot(M2.H4_time./86400, M2.H4_wh_p.*Pa2MPa, 'k-', 'LineWidth', 1);
        xlabel('days'); ylabel('well head pressure (MPa)')
        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'k';
        drawnow
    end
end

if saveflag 
    save('injection_output', 'r', 't_s', 'p_s', 'pWH', 'q_s');     
end


end