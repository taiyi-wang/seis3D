% Benchmark linear, radial pressure diffusion solver against analytical
% solution and global mass balance. This benchmark neglects well
% storativity and viscous pressure loss.

%% add paths
current_dir  = pwd;
idcs   = strfind(current_dir ,'/');
above_dir = current_dir(1:idcs(end)-1); % one level above current directory
addpath(fullfile(above_dir, 'source_code/.'));
%% Set parameters
testflag = 1;    % 1 for testing
rglgridflag = 0; % 0 for using log grid, 1 for using regular grid
plotflag = 0;    % 0 for not plotting results for each injection experiment
saveflag = 0;    % 0 for not saving results into .mat files

% unit conversion
Pa2MPa = 1e-6;

% fault zone properties
M.w = 4; % width (m)
M.k = 1e-12; % permeability (m^2)
M.beta_phi = 0.05*1e-8; % porosity*(fluid+pore compressibilities) (1/Pa)
M.eta_v = 1e-3; % fluid viscosity (Pa*s)
M.p_pore = 0;   

% well properties
M.R_w = 0.1778/2; % radius (m)
M.H_w = 0;        % well height (m)
M.beta_w = 0;     % well compressibility (Pa^-1)
M.f = 0; % Darcy-Weisbach friction factor
M.rho = 1000; % fluid density (kg/m^3)
M.g = 9.8;

% radial distance, solve for rmin<r<rmax
M.rmin = M.R_w; % minimum radius (m)
M.rmax = 10e3; % maximum radius (m)
M.nr = 100; % number of grid intervals = number of grid points - 1

% simulation set up
M.duration = 300;
M.Q0 = 0.1; % step function injection rate (m^3/s)
M.dt_max = M.duration./10;

% refinement factor (to easily refine in both space and time)
refine_list = 1:5;
r_store = cell(length(refine_list), 1);
t_store = cell(length(refine_list), 1);
q_store = cell(length(refine_list), 1);
p_store = cell(length(refine_list), 1);

%% Run simulation at various spatial/temporal refinement levels
nr = M.nr;
dt_max = M.dt_max;
for i = 1:length(refine_list)
    refine = refine_list(i);
    M.nr = nr*refine;
    M.dt_max = dt_max/refine;
    [r, t_s, p_s, ~, ~, q_s] = injection2(M, testflag, rglgridflag, plotflag, saveflag);

    r_store{i} = r;
    t_store{i} = t_s;
    q_store{i} = q_s;
    p_store{i} = p_s;
    
end
%% convergence of pressure time series with mesh refinement
% 1. convergence at injector (well bottom) --------------------------------
alpha_list = linspace(0.1, 1, length(refine_list)); % a list of different levels of transparency for plotting

figure;
for i = 1:length(refine_list)
    t = t_store{i};
    P = p_store{i}; % pressure history at all locations
    p = P(1, :);
    plt = plot(t./3600, p.*Pa2MPa, 'k-', 'LineWidth', 1); 
    plt.Color(4) = alpha_list(i); hold on;
    xlabel('hour (s)'); ylabel('pressure change (MPa)')
end
title('r = R_w')

% 2. convergence away from the injector at 10 x well radius ---------------
% extract the highest resolution results
r_fine = r_store{end};
p_fine = p_store{end};
t_fine = t_store{end};

% analytical solution 
NRw = 1e1; % number of well radii away from the center of injection
r0 = M.R_w*NRw; 
phi = 0.05; 
beta = M.beta_phi/phi;
p_t_alyt = alyt_step_p_sol(t_fine, r0, M.Q0, M.k, beta, phi, M.eta_v, M.w, 0);

% interpolate numerical solution at r0 at different refinement levels
p_interp_levels = cell(length(refine_list), 1);
for i = 1:length(refine_list)
    t_list = t_store{i};
    for j = 1:size(t_list, 2)
        p_array = p_store{i};
        p_interp(j) = interp1(r_store{i}, p_array(:, j), r0);
    end
    p_interp_levels{i} = p_interp;
end

figure;
for i = 1:length(refine_list)
    plt = plot(t_store{i}./3600, p_interp_levels{i}.*Pa2MPa, 'k-', 'LineWidth', 1); 
    plt.Color(4) = alpha_list(i); hold on;
end
plot(t_store{end}./3600, p_t_alyt.*Pa2MPa, 'r-', 'LineWidth', 1); 
xlabel('hours'); ylabel('pressure change (MPa)'); legend('numerical', 'analytical')
title(sprintf('r = %i x R_w', NRw))

% 3. pressure profile over time ----------------------------------------------
figure;
r = r_store{end};
P = p_store{end}; % pressure history at all locations
plot(r,P.*Pa2MPa, 'k-', 'LineWidth', 1); 
xlim([0, 400]);
xlabel('radial distance from injector (m)'); ylabel('pressure change (MPa)')
title('domain R_{max} = 10km')

%% global conservation of mass
% run simulation using no flow boundary condition
% spatial discretization
refine = 20; %max(refine_list);
nr = 100; % number of grid intervals = number of grid points - 1

M.nr = nr*refine;
M.dt_max = dt_max/refine;
[r, t_s, p_s, pWH, dpw0, q_s] = injection1(M, testflag, rglgridflag, plotflag, saveflag);

% compute numerical average pressure at every time step
p_avg = nan(length(t_s), 1);
R = max(r) - min(r);
xub = 5e2; xlb = -5e2;
yub = 5e2; ylb = -5e2;
x_rg = linspace(xlb, xub, 100);
y_rg = linspace(ylb, yub, 100);
[X_rg, Y_rg] = meshgrid(x_rg, y_rg); % regular Cartesian grid inset into the radial domain
R_rg = sqrt(X_rg.^2 + Y_rg.^2);       % radial distance corresponding to regular grid
A = (xub - xlb)*(yub - ylb);
for i = 1:length(t_s)
    % interpolate the p(r, theta) onto regular, Cartesian grid p(x, y)
    for m = 1:length(x_rg)
        for n = 1:length(y_rg)
            p_xy(n, m) = interp1(r, p_s(:, i), R_rg(n, m));
        end
    end
    
    % compute spatially averaged pressure over the rectangular domain
    p_avg(i) = trapz(x_rg, trapz(y_rg, p_xy, 1),2)/A; % spatially averaged pressure
end

% analytical solution
p_alyt = M.Q0.*t_s./(M.w*M.beta_phi*A);

figure;
yyaxis left;
plot(t_s./3600, p_avg.*Pa2MPa, 'k--', 'LineWidth', 2); hold on;
plot(t_s./3600, p_alyt.*Pa2MPa, 'r-', 'LineWidth', 2);
ylabel('pressure change (MPa)');
yyaxis right;
plot(t_s./3600, q_s, '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
ylabel('injection rate (m^3/s)'); ylim([0.05, M.Q0+0.05]);
xlabel('time (hours)'); 
legend('numerical', 'analytical', 'injection input')
title(sprintf('pressure averaged over %i x %i m domain', (xub - xlb), (yub - ylb)));
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%% storativity
M.k = 0;                % non-permeable fault
M.H_w = 1000;           % well height
M.beta_w = 4.41e-10;    % well compressibility
%{
[r, t_s, p_s, pWH, t_s_w0, dpw0, q_s] = injection2(M, testflag, rglgridflag, plotflag, saveflag);

figure;
plot(t_s_w0./3600, dpw0); hold on;
plot(t_s./3600, p_s(1,:)-p_s(1, 1));
legend('reference no flow solution', 'zero-permeability solution')
xlabel('hours')
%}

[r, t_s, p_s, pWH, dpw0, q_s] = injection1(M, testflag, rglgridflag, plotflag, saveflag);
figure;
plot(t_s./3600, dpw0); hold on;
plot(t_s./3600, p_s(1,:)-p_s(1, 1));
legend('reference no flow solution', 'zero-permeability solution')
xlabel('hours')







