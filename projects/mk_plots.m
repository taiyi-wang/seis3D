% Make plots for publication
% Note: 
% 1. the code is meant to be run sequentially

% Taiyi Wang 11/30/2021

current_dir  = pwd;
idcs   = strfind(current_dir ,'/');
above_dir = current_dir(1:idcs(end)-1); % one level above current directory

%addpath(fullfile(above_dir, '/output'))
addpath(fullfile(above_dir ,'/source_code'))
addpath(fullfile(above_dir ,'/output/final/.'))
addpath(fullfile(above_dir ,'/input'))

%% Adjust font and line
set(0,'DefaultTextFontSize',16)
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultLineLineWidth', 2);
set(0,'DefaultLineMarkerSize',12);

%% Make schematic to show the main fault and seismic patches
load('M1.mat'); 
load('M2.mat');

% unpack variables
ss_locs = M1.ss_locs;
ss_dx = M1.ss_dx;
ss_dz = M1.ss_dz;
fault_depth = M1.fault_depth;

% define rectangle's vertices locations counterclockwise from upper left
ss_X = [ss_locs(:, 1)'-ss_dx/2; ss_locs(:, 1)'-ss_dx/2; ss_locs(:, 1)'+ss_dx/2; ss_locs(:, 1)'+ss_dx/2];
ss_Z = [ss_locs(:, 2)'+ss_dz/2; ss_locs(:, 2)'-ss_dz/2; ss_locs(:, 2)'-ss_dz/2; ss_locs(:, 2)'+ss_dz/2];
ss_Y = repmat(-ss_locs(:, 3)', [4, 1]);

% define main fault's verticies (not reflective of the true dimension of main fault)
as_X = [-1.5e3; -1.5e3; 1.5e3; 1.5e3];
as_Z = [1.5e3; -1.5e3; -1.5e3; 1.5e3];
as_Y = [-fault_depth; -fault_depth; -fault_depth; -fault_depth];

% 1. plot locations of seismic patches/main fault and well
figure;
plt1 = fill3(ss_X./1e3, ss_Z./1e3, ss_Y./1e3, [1, 0.4, 0]); hold on;
set(plt1,'edgecolor', [1, 0.4, 0]);
fill3(as_X./1e3, as_Z./1e3, as_Y./1e3, 'k');
plot3(repelem(0, 5, 1), repelem(0, 5, 1), linspace(-4.1, -4.085, 5), 'Color', [0.2, 0.4, 0.8], 'LineWidth', 2);
grid on;
alpha(0.3);
xlabel('x, along dip (km)'); ylabel('z, along strike (km)'); zlabel('y, depth (km)')
xlim([-1, 1]); ylim([-1, 1]); 

% 2. plot histogram of depth distribution
figure;
histogram(ss_locs(:, 3)./1000, 100, 'FaceColor', [1, 0.4, 0], 'FaceAlpha', 0.8, 'LineWidth', 1);
xlabel('depth (km)')


%% Plot injection rate, well head pressure data, with prediction
load('injection_output.mat');

% unpack variables
x0_4 = M2.x0_4; z0_4 = M2.z0_4; X = M1.X; Z = M1.Z;
H4_time = M2.H4_time; H4_wh_p = M2.H4_wh_p;

% observation point coordinates
xx_0_idx = find(X(1,:) == x0_4);
zz_0_idx = find(Z(:,1) == z0_4);
wh_time = H4_time;
p_wh_m = H4_wh_p./1e6;

figure;
yyaxis left
plot(t_s./86400, q_s.*1000, '-', 'Color', [0.5, 0.5, 0.5]);
ylabel('Q (L/s)')
yyaxis right
plot(wh_time./86400, p_wh_m-p_wh_m(1), 'Color', 'k'); hold on;
plot(t_s./86400, (pWH-pWH(1))./1e6, 'Color', [0 0.4470 0.7410], 'LineStyle', '-')
ylabel('\delta p_{wh} (MPa)'); ylim([-1, 30]);
xlabel('days since injection started')
ax = gca;
ax.YAxis(1).Color = [0.5, 0.5, 0.5];
ax.YAxis(2).Color = [0 0.4470 0.7410];

%% Pore pressure vs. normal stress at injector
figure;
plot(t_s./86400, repelem(M1.tauy_as_0(1,1), length(t_s))./1e6, 'k-'); hold on;
plot(t_s./86400, p_s(:, 1)./1e6, '-');
xlabel('days'); ylabel('MPa');
legend('lithostatic + far-field', 'fault pore pressure at injector')

%% Compare shear stress vs. pressure perturbation evolution
load('as_output.mat', 't_as_s', 'p_as_s', 'taux_as_s')
X = M1.X; Z = M1.Z;

Nt = length(t_as_s);                                % number of time steps
x_locs = X(1,:); z_locs = Z(:,1);

t_1day_interval = 1:2:(t_as_s(end)/86400); % equally spaced time (days)
t_days = t_as_s./86400;                  % time from simulation (days)

% interpolate dp, dtaux, for evenly spaced time
new_t = t_as_s(1):86400:t_as_s(end); % time with one day intervals
new_dp_z = zeros(length(z_locs), length(new_t)); % pressure change along z
new_dp_x = zeros(length(x_locs), length(new_t)); % pressure change along x
new_dtaux_z = zeros(length(z_locs), length(new_t)); % shear stress change along z
new_dtaux_x = zeros(length(x_locs), length(new_t)); % shear stress change along x

% a profile along z, passing through the well
for i = 1:length(z_locs)
    ts_z = squeeze(p_as_s(i, xx_0_idx, :)); % time series at a specific point along z
    new_dp = interp1(t_as_s, ts_z, new_t) - ts_z(1);
    new_dp_z(i, :) = new_dp;
    
    ts_z = squeeze(taux_as_s(i, xx_0_idx, :)); % time series at a specific point along z
    new_dtaux = interp1(t_as_s, ts_z, new_t) - ts_z(1);
    new_dtaux_z(i, :) = new_dtaux;
end

% a profile along x, passing through the well
for i = 1:length(x_locs)
    ts_x = squeeze(p_as_s(zz_0_idx, i, :)); % time series at a specific point along x
    new_dp = interp1(t_as_s, ts_x, new_t) - ts_x(1);
    new_dp_x(i, :) = new_dp;
    
    ts_x = squeeze(taux_as_s(zz_0_idx, i, :)); % time series at a specific point along x
    new_dtaux = interp1(t_as_s, ts_x, new_t) - ts_x(1);
    new_dtaux_x(i, :) = new_dtaux;
end

clear p_as_s taux_as_s

%%
figure;
for i = 1:2:length(new_t)
    plot(Z(:, xx_0_idx)'./1000, new_dp_x(:,i)./1e6, 'r', 'LineWidth', 1);hold on;  
    plot(X(zz_0_idx,:)./1000, new_dp_z(:,i)./1e6, 'r--', 'LineWidth', 1)
end
xline(110/1e3, 'k--', 'LineWidth', 1); xline(-110/1e3, 'k--', 'LineWidth', 1);
xlim([-3, 3]);
xlabel('locations (km)'); ylabel('pressure change (MPa)');
legend('along x-axis (along dip)', 'along z-axis (along strike)');
title('time interval = 2 day');

figure;
for i = 1:2:length(new_t)
    plot(Z(:, xx_0_idx)'./1000, new_dtaux_x(:,i)./1e6, 'b', 'LineWidth', 1);hold on;  
    plot(X(zz_0_idx,:)./1000, new_dtaux_z(:,i)./1e6, 'b--', 'LineWidth', 1)
end
xlim([-3, 3]);
xlabel('locations (km)'); ylabel('shear stress change (MPa)');
legend('along x-axis (along dip)', 'along z-axis (along strike)');
title('time interval = 2 day');

%% Plot cumulative aseismic slip along each fault dimension
load('as_output.mat', 'Dx_as_s')

x_locs = X(1,:);
z_locs = Z(:,1);

% first interpolate displacement for evenly spaced time
new_t = t_as_s(1):86400:t_as_s(end); % time with one day intervals
new_D_z = zeros(length(z_locs), length(new_t)); % displacement along z
new_D_x = zeros(length(x_locs), length(new_t)); % displacement along x

% a profile along z, passing through the well
for i = 1:length(z_locs)
    ts_z = squeeze(Dx_as_s(i, xx_0_idx, :)); % time series at a specific point along z
    new_d = interp1(t_as_s, ts_z, new_t);
    new_D_z(i, :) = new_d;
end

% a profile along x, passing through the well
for i = 1:length(x_locs)
    ts_x = squeeze(Dx_as_s(zz_0_idx, i, :)); % time series at a specific point along x
    new_d = interp1(t_as_s, ts_x, new_t);
    new_D_x(i, :) = new_d;
end

figure;
for i = 1:2:length(new_t)
    plot(Z(:, xx_0_idx)'./1000, new_D_x(:,i).*100, 'k', 'LineWidth', 1);hold on;  
    plot(X(zz_0_idx,:)./1000, new_D_z(:,i).*100, 'k--', 'LineWidth', 1)
end
xlim([-3, 3]);
xlabel('locations (km)'); ylabel('aseismic slip (cm)');
legend('along x-axis (along dip)', 'along z-axis (along strike)');
title('time interval = 2 day');

clear Dx_as_s

%% Plot seismic slip along each dimension
load('ss_output_15.mat', 'Dx_ss_s')

% find indices of spring sliders within thin strips along x and z passing
% through injection point 
ss_x_idc = find(abs(ss_locs(:,2)) <= 100); % along x
ss_z_idc = find(abs(ss_locs(:,1)) <= 100); % along z

ss_xs = ss_locs(ss_x_idc,1);
ss_zs = ss_locs(ss_z_idc,2);

for i = 1:length(ss_x_idc)
    ss_x_idx = ss_x_idc(i);
    ss_D_x(i) = Dx_ss_s{ss_x_idx}(end);
end

for i = 1:length(ss_z_idc)
    ss_z_idx = ss_z_idc(i);
    ss_D_z(i) = Dx_ss_s{ss_z_idx}(end);
end

figure;
plot(ss_xs./1000, ss_D_x.*100 , 'k.', 'MarkerSize', 30); hold on;
plot(ss_zs./1000, ss_D_z.*100 , 'ko', 'MarkerSize', 8); 
xlim([-1, 1]); legend('along x-axis (along dip)', 'along z-axis (along strike)')
xlabel('location (km)'); ylabel('final seismic slip (cm)');

%% Plot H4 aseismic slip and (one way coupled) seismic events
load('as_output.mat', 'Vx_as_s')
load('ss_output_15.mat')
% Unpack variables
Nx = M1.Nx; Nz = M1.Nz; p_pore = M2.p_pore;
N_ss = M1.N_ss; mu = M1.mu;

mid_x_idx = (Nx+1)/2;
mid_z_idx = (Nz+1)/2;

% Velocity on main fault in log10(m/s)
Vx_as = squeeze(Vx_as_s(:,mid_x_idx,:));
Vx_as_max = max(Vx_as, [],'all');
Vx_as_min = min(Vx_as, [],'all');

z_locs = Z(:, 1);
xlocs_ss = ss_locs(:,1)./1000; % x coordinates of the spring sliders in km
zlocs_ss = ss_locs(:,2)./1000; % x coordinates of the spring sliders in km
[T_grid, Z_grid] = meshgrid(t_as_s, x_locs);

% retrieve locations of observed seismic events along z = 0 km, for x >= 0 km
load('Cooper_Basin_Catalog_HAB_4.mat')
event_lat = Catalog(4).val; event_lon = Catalog(5).val; elevation = Catalog(6).val; event_t = Catalog(2).val;
event_t = (event_t-event_t(1)); event_Mw = Catalog(8).val; event_M0 = Catalog(7).val;
HB4_lat = -27.8115; HB4_lon = 140.7596;% from well completion report

% convert event locations to meters with regard to HB4
event_local = llh2local([event_lon'; event_lat'; elevation'],[HB4_lon; HB4_lat; 0]);
event_x = event_local(1,:).*1000; event_z = event_local(2,:).*1000; 

% find events along x-axis (west-east), extending from HB4 towards west-east
event_narroaw_idx = intersect(find(event_x<100 & event_x>-100), find(event_Mw > 0));

figure;
p = pcolor(T_grid./86400, Z_grid./1000, log10(Vx_as)); colormap(parula);colorbar; caxis([log10(Vx_as_min), log10(Vx_as_max)-2]);
p.FaceColor = 'interp'; set(p, 'EdgeColor', 'none'); p.FaceAlpha = 1; hold on;

for i = 1:length(event_narroaw_idx)
    plt = plot(event_t(event_narroaw_idx(i))-3, event_z(event_narroaw_idx(i))./1000, 'w.', 'MarkerSize', event_Mw(event_narroaw_idx(i))*10);% seismic catalogue starts on 11/10/2012, whereas pressure data starts on 11/13/2012
    plt.Color(4) = 0.2; 
end

quakes = find_quakes(t_ss_s, Vx_ss_s, Dx_ss_s, M1);

for i = 1:size(quakes, 1)
    plot(quakes{i, 2}, zlocs_ss(i), '.', 'Color', [0.8500, 0.3250, 0.0980, 0.2], 'MarkerSize', abs(quakes{i, 4}).*10);  
end

% make a legend for seismic event moment magnitude
%plot(1.1, -0.6, 'w.', 'MarkerSize', 2*10); %text(1.2, -0.5, 'M_w = 2', 'Color', 'w');
%plot(1.1, -0.75, 'w.', 'MarkerSize', 1.5*10); %text(1.2, -0.65, 'M_w = 1.5', 'Color', 'w');
%plot(1.1, -0.9, 'w.', 'MarkerSize', 1*10); %text(1.2, -0.8, 'M_w = 1', 'Color', 'w');
view([0, 90]);
xlim([min(t_s./86400), 17]); ylim([-1, 1]);
ylabel('z (km)'); xlabel('time (days)'); %title('Color scale = aseismic slip velocity (cm/day); contour = perturbation pressure (MPa)')

clear Vx_as_s

%% Histograms of seismicity
% 1. time dependence of seismicity
data_event_t = event_t-3; % seismic catalogue starts on 11/10/2012, whereas pressure data starts on 11/13/2012
data_event_t = data_event_t(data_event_t<=17);
ss_event_t = cell2mat(quakes(:,2));

figure;
histogram(data_event_t, 'Normalization', 'Probability', 'FaceAlpha', 0.5, 'FaceColor', [1, 1, 1], 'LineWidth', 1, 'BinWidth', 0.2); hold on;
histogram(ss_event_t, 'Normalization', 'Probability', 'LineWidth', 1, 'FaceAlpha', 0.5, 'BinWidth', 0.2);
%histogram(data_event_t, 'Normalization', 'cdf', 'FaceColor', 'none', 'LineWidth', 1, 'BinWidth', 0.2, 'DisplayStyle', 'stairs'); hold on;
%histogram(ss_event_t, 'Normalization', 'cdf', 'LineWidth', 1, 'BinWidth', 0.2, 'DisplayStyle', 'stairs');

%[f_d,xi_d] = ksdensity(data_event_t, 'NumPoints',1000);
%[f_p,xi_p] = ksdensity(ss_event_t, 'NumPoints',1000);
%plot(xi_d, f_d, 'Color', 'k', 'LineWidth', 1);
%plot(xi_p, f_p, 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1);
%legend('catalogue', 'simulated')
%xlabel('days since injection started'); ylabel('frequency of seismic events');
xlim([0, 18]);

% 2. seismicity as a function of spatial coordinates

% first find the x and z coordinates of predicted earthquakes
pred_x = [];
pred_z = [];

for i = 1:size(quakes,1)
    Nevent = length(quakes{i, 2});
    pred_x = [pred_x;ss_locs(repelem(quakes{i, 1}, Nevent), 1)];
    pred_z = [pred_z;ss_locs(repelem(quakes{i, 1}, Nevent), 2)];
end

% find the indices of earthquakes within narrow bands along each direction
data_narroaw_z_idx = intersect(find(event_z<100 & event_z>-100), find(event_Mw > 0));
data_narroaw_x_idx = intersect(find(event_x<100 & event_x>-100), find(event_Mw > 0));

pred_narroaw_z_idx = intersect(find(pred_z<500 & pred_z>-500), find(event_Mw > 0));
pred_narroaw_x_idx = intersect(find(pred_x<500 & pred_x>-500), find(event_Mw > 0));

figure;
histogram(event_x(data_narroaw_z_idx)./1000, 'Normalization', 'Probability', 'FaceAlpha', 0.5, 'FaceColor', [1, 1, 1], 'LineWidth', 1, 'BinWidth', 0.03);hold on;
histogram(pred_x(pred_narroaw_z_idx)./1000, 30, 'Normalization', 'Probability', 'LineWidth', 1, 'FaceAlpha', 0.5, 'BinWidth', 0.03); 
xlim([-1.5, 1.5]);
xlabel('distance along dip (km)');
ylabel('normalized frequency')

figure;
histogram(event_z(data_narroaw_x_idx)./1000, 'Normalization', 'Probability', 'FaceAlpha', 0.5, 'FaceColor', [1, 1, 1], 'LineWidth', 1, 'BinWidth', 0.03);hold on;
histogram(pred_z(pred_narroaw_x_idx)./1000, 30, 'Normalization', 'Probability', 'LineWidth', 1, 'FaceAlpha', 0.5, 'BinWidth', 0.03); 
xlim([-1.5, 1.5]);
xlabel('distance along strike (km)');
ylabel('normalized frequency')

%% migration of seismicity over time
quake_x = ss_locs(cell2mat(quakes(:, 1)), 1);
quake_z = ss_locs(cell2mat(quakes(:, 1)), 2);
quake_r = sqrt(quake_x.^2 + quake_z.^2);

event_r = sqrt(event_x.^2 + event_z.^2);

figure;
scatter(event_t-3, event_r./1e3, 50, 'k', 'filled', 'MarkerFaceAlpha', 0.3); hold on;

for i = 1:length(quake_r)
    quake_ts = quakes{i, 2};
    for j = 1:length(quake_ts)
        quake_t = quake_ts(j);
        scatter(quake_t, quake_r(i)./1e3, 50, 'r', 'filled'); hold on;
    end
end
xlim([0, 17]);


xlabel('days since injection began'); ylabel('distance from injector (km)')
clear quake_x quake_z quake_r

%% Normalized cumulative seismic/aseismic moment release
load('as_output.mat', 'Dx_as_s')

A = M1.dx .* M1.dz; % area of a grid element on main fault
A_ss = M1.ss_dx*M1.ss_dz;
mu = M1.mu;

% 1. compare data/prediction in seism ic moment release
t_m = (0:0.05:17).*86400; % evenly spaced time for resampling

% compute predicted total seismic moment release (default case)
cM_ss_p = zeros(1, length(t_m));
load('ss_output_15.mat', 't_ss_s', 'Dx_ss_s')
for i = 1:N_ss
    t = t_ss_s{i};
    d = Dx_ss_s{i};
    m = A_ss .* mu .*  d; % cumulative moment
    
    cM_ss_p = cM_ss_p + interp1(t, m, t_m);
end
clear t d m t_ss_s Dx_ss_s

% compute predicted total seismic moment release (12 MPa initial shear stress)
cM_ss_p_12 = zeros(1, length(t_m));
load('ss_output_12.mat', 't_ss_s', 'Dx_ss_s')
for i = 1:N_ss
    t = t_ss_s{i};
    d = Dx_ss_s{i};
    m = A_ss .* mu .*  d; % cumulative moment
    
    cM_ss_p_12 = cM_ss_p_12 + interp1(t, m, t_m);
end
clear t d m t_ss_s Dx_ss_s

% compute predicted total seismic moment release (17 MPa initial shear stress)
cM_ss_p_17 = zeros(1, length(t_m));
load('ss_output_17.mat', 't_ss_s', 'Dx_ss_s')
for i = 1:N_ss
    t = t_ss_s{i};
    d = Dx_ss_s{i};
    m = A_ss .* mu .*  d; % cumulative moment
    
    cM_ss_p_17 = cM_ss_p_17 + interp1(t, m, t_m);
end
clear t d m t_ss_s Dx_ss_s

% compute observed total seismic moment release
cM_ss_d = cumsum(event_M0);
[~, idx_d] = min(abs(event_t - 3 - 17));

figure;
plot(event_t-3, cM_ss_d./cM_ss_d(idx_d), 'k-'); hold on;
plot(t_m./86400, cM_ss_p./max(cM_ss_p), '-');
plot(t_m./86400, cM_ss_p_12./max(cM_ss_p_12), '--');
plot(t_m./86400, cM_ss_p_17./max(cM_ss_p_17), ':');
xlim([0, 17]);
xlabel('days'); ylabel('normalized seismic cumulative moment')
legend('catalogue', '\tau_{yx}^0 = 15 MPa', '\tau_{yx}^0 = 12 MPa', '\tau_{yx}^0 = 17 MPa')

% 2. compare predicted aseismic vs. observed seismic moment release

% compute predicted total aseismic moment release
cM_as_p = squeeze(sum(sum(Dx_as_s .* mu .* A, 1), 2));

figure;
plot(t_as_s./86400, cM_as_p, 'k-'); hold on;
plot(event_t-3, cM_ss_d, 'k--', 'MarkerSize', 1); 
xlim([0, 17]); xlabel('days'); ylabel('normalized cumulative moment')
legend('simulated cumulative aseismic moment', 'observed cumulative seismic moment')

clear idx_d

clear Dx_as_s

%% Main fault phase diagram and velocity time series
load('as_output.mat', 'Vx_as_s', 'psi_as_s')
% pick a location 
xidx = 108; zidx = 108;
a_as = M1.a(zidx, xidx); b_as = M1.b(zidx, xidx);
f0 = M1.f0; V0 = M1.V0;
Vx_as = squeeze(Vx_as_s(zidx, xidx, :));
psi_as = squeeze(psi_as_s(zidx, xidx, :));

% steady state 
f_ss = @(v) f0 + (a_as-b_as) * log(v/V0);

% frictional coefficient
f = a_as .* asinh(Vx_as/(2*V0).*exp(psi_as/a_as));

figure;
plot(log(Vx_as)- log(V0), f, 'k', 'LineWidth',1); hold on;  % phase diagram
plot(log(Vx_as) - log(V0), f_ss(Vx_as), 'r-', 'LineWidth', 1);  % steady state
plot(log(Vx_as) - log(V0), a_as.*log(Vx_as)+0.96 , 'b-', 'LineWidth', 1);
plot(log(Vx_as(1)) - log(V0), f(1), 'k.', 'MarkerSize', 20); % starting point
xlabel('log_{10} v_x/v_0'), ylabel('f');

figure;
semilogy(t_as_s, Vx_as, 'k-'); 
xlabel('time (s)'); ylabel('log_{10} v');

clear Vx_as_s psi_as_s
%% spring slider phase diagram
load('ss_output_15.mat', '')
% Note: spring slider # 86 has 2 events
% choose the index of spring slider
choose_idc = 18;

% unpack variables
a_ss = M1.a_ss; b_ss = M1.b_ss;
f0 = M1.f0; V0 = M1.V0;

% steady state 
f_ss = @(v) f0 + (a_ss-b_ss) * log(v/V0);

% plot selected spring sliders on the same phase diagram
figure;
for i = 1:length(choose_idc)
    choose_idx = choose_idc(i);
    
    % compute the change in frictional coefficient
    f = a_ss .* asinh(Vx_ss_s{choose_idx}/(2*V0).*exp(Psi_ss_s{choose_idx}/a_ss));
    
    semilogx(Vx_ss_s{choose_idx}, f, 'k', 'LineWidth',1);hold on;  % phase diagram
    semilogx(Vx_ss_s{choose_idx}(1), f(1), 'k.', 'MarkerSize', 20); % starting point
end
semilogx(10.^linspace(-15, 1, 100), f_ss(10.^linspace(-15, 1, 100)), 'r', 'LineWidth', 1);   % steady state frictional coefficient
xlabel('log_{10} v_x'), ylabel('f = \tau_{yx}/(\tau_{yy} - p_{pore})'); hold off;
clear choose_idx

figure;
for i = 1:length(choose_idc)
    choose_idx = choose_idc(i);
    
    semilogy(t_ss_s{choose_idx}./86400, Vx_ss_s{choose_idx}, 'k-'); hold on;
end

%% spring slider pressure and aseismic loading stress history
choose_idx = 18;
% compute total strength change
f = a_ss .* asinh(Vx_ss_s{choose_idx}/(2*V0).*exp(Psi_ss_s{choose_idx}/a_ss));
f0 = f(1);
qd_stress_change = taux_ss_s{choose_idx} - taux_ss_s{choose_idx}(1); % quasi-dynamic stress change
strength_change = f.*qd_stress_change;

% compute stress change due to aseismic loading
shear_stress_change = dtaux_ss(:, choose_idx);
normal_stress_change = dtauy_ss(:, choose_idx);

figure;
subplot(2, 1, 1)
semilogy(t_ss_s{choose_idx}./86400, Vx_ss_s{choose_idx}, 'k-'); 
ylabel('log_{10} v_x')
subplot(2, 1, 2)
plot(t_ss_s{choose_idx}./86400, strength_change./1e6, 'k'); ylabel('(MPa)'); hold on; % total shear strength change
plot(t_as_s./86400, shear_stress_change(2:end)./1e6, 'b--'); % shear stress change due to aseismic loading
plot(t_as_s./86400, f0.*normal_stress_change(2:end)./1e6, 'b'); % shear strength change due to normal stress change produced by aseismic loading
plot(t_ss_s{choose_idx}./86400, -f0.*(p_ss_s{choose_idx} - p_ss_s{choose_idx}(1))./1e6, 'r'); % shear strength change due to normal stress change produced by pressure diffusion
legend('\delta \tau_{str}', '\delta \tau_{as}', 'f_0 \cdot \delta \sigma_{as}', '-f_0 \cdot \delta p_{ss}')
xlabel('time (days)')

%% Plot observed seismic catalogue
% retrieve locations of observed seismic events along z = 0 km, for x >= 0 km
load('Cooper_Basin_Catalog_HAB_4.mat')
event_lat = Catalog(4).val; event_lon = Catalog(5).val; elevation = Catalog(6).val; event_t = Catalog(2).val;
event_t = (event_t-event_t(1)); event_Mw = Catalog(8).val; event_M0 = Catalog(7).val;
HB4_lat = -27.8115; HB4_lon = 140.7596;% from well completion report

% convert event locations to meters with regard to HB4
event_local = llh2local([event_lon'; event_lat'; elevation'],[HB4_lon; HB4_lat; 0]);
event_x = event_local(1,:).*1000; event_z = event_local(2,:).*1000; 

figure;
scatter3(event_x./1e3, event_z./1e3, elevation', (event_Mw./max(event_Mw)+0.1).*5e2, event_t, '.');hold on; colorbar; 
plot3(repelem(0, 5, 1), repelem(0, 5, 1), linspace(-3.5, -4.085, 5), 'Color', [0.2, 0.4, 0.8], 'LineWidth', 2);
xlabel('East (km)'); ylabel('North (km)');
colormap('jet'); axis equal;
grid on;

%% Plot histogram of strength drop vs. shear stress change due to aseismic loading
[Nseismic, ~] = size(quakes); % number of patches that became seismic
count = 0;
for i = 1:Nseismic
    idx = quakes{i, 1}; % index of seismic patches
    tes = quakes{i, 2}.*86400; % timing of earthquakes (s)
    ps = p_ss_s{idx};   % pressure history at the seismic patches
    ts = t_ss_s{idx};   % time of each seismic patch
    p0 = ps(1);         % initial pressure at this seismic patch
    dtauxs = dtaux_ss(:, idx);    % shear stress change due to aseismic loading
    dtauys = dtauy_ss(:, idx);    % normal stress change due to aseismic loading
    f = M1.a_ss .* asinh(Vx_ss_s{idx}/(2*M1.V0).*exp(Psi_ss_s{idx}/M1.a_ss));
    f0 = f(1);
    for j = 1:length(tes)
        count = count + 1;
        te = tes(j);
        str_drop = interp1(ts, -(ps - p0).*f0, te);
        dtau_x = interp1(t_as_s, dtauxs(2:end), te);
        dtau_y = interp1(t_as_s, dtauys(2:end), te);

        str_drops(count) = str_drop;
        dtau_xs(count) = dtau_x;
        f0dtau_ys(count) = dtau_y.*f0;
        quakex(count) = ss_locs(idx, 1);
        quakez(count) = ss_locs(idx, 2);
        quaket(count) = te;
    end
end

figure;
histogram(dtau_xs./1e6, 'BinWidth', 0.1, 'FaceColor', 'b', 'EdgeAlpha', 1, 'FaceAlpha', 0.5, 'LineWidth', 1);hold on;
histogram(str_drops./1e6, 'BinWidth', 0.1,'FaceColor', 'r', 'EdgeAlpha', 1, 'FaceAlpha', 0.5, 'LineWidth', 1); 
histogram(f0dtau_ys./1e6, 'BinWidth', 0.1, 'FaceColor', [0.9290 0.6940 0.1250], 'EdgeAlpha', 1,'FaceAlpha', 0.5, 'LineWidth', 1); 
%legend('\Delta \tau_{as}', '- f_0 \cdot aiyi\Delta p ', 'f_0 \cdot \Delta \sigma_{as}')
xlabel('MPa'); ylabel('count'); xlim([-4, 4]); 

% reorder pressure change according to distance from injector
[~, sortidx] = sort(sqrt(quakex.^2 + quakez.^2));
quaker = sqrt(quakex.^2 + quakez.^2);
Nevents = length(sortidx);
%figure
scatter(dtau_xs(sortidx)./1e6, (1:Nevents)./2.3, 100, 'b', 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'k'); hold on;
scatter(str_drops(sortidx)./1e6, (1:Nevents)./2.3, 100, 'r', 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'k');
scatter(f0dtau_ys(sortidx)./1e6, (1:Nevents)./2.3, 100, [0.9290 0.6940 0.1250], 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'k');
ax=gca;ax.LineWidth=1; 
ylim([0, 100]);

% plot the location of aseismic loading triggered versus fault weakening
% triggered events
as_trigger_idc = union(find(abs(dtau_xs) >= abs(str_drops)), find((abs(f0dtau_ys) >= abs(str_drops))));
wk_trigger_idc = intersect(find(abs(str_drops)>=abs(dtau_xs)), find(abs(str_drops)>=abs(f0dtau_ys)));

[Xq, Zq] = meshgrid(linspace(min(X(:)), max(X(:)), 5e3), linspace(min(Z(:)), max(Z(:)), 5e3));
Dx_end = Dx_as_s(:,:,end);
Dx_end = interp2(X, Z, Dx_end, Xq, Zq);

figure;
contour(Xq./1000, Zq./1000, Dx_end, 'LineWidth', 2, 'Color', 'k', 'ShowText','on'); hold on;
plot(quakex(wk_trigger_idc)./1e3, quakez(wk_trigger_idc)./1e3, 'r.', 'MarkerSize', 20);
plot(quakex(as_trigger_idc)./1e3, quakez(as_trigger_idc)./1e3, 'b+', 'MarkerSize', 10);
plot(ss_locs(choose_idx, 1)./1e3, ss_locs(choose_idx, 2)./1e3, 'r+', 'MarkerSize', 20); 
axis equal
xlim([-1.5, 1.5]); ylim([-1.5, 1.5]);

%% plot space-time of triggering mechanism
quaker_wk = sqrt(quakex(wk_trigger_idc).^2 + quakez(wk_trigger_idc).^2);
quaker_as = sqrt(quakex(as_trigger_idc).^2 + quakez(as_trigger_idc).^2);
quaket_wk = quaket(wk_trigger_idc);
quaket_as = quaket(as_trigger_idc);

% radial distance of chosen spring slider
chosen_r = sqrt(ss_locs(choose_idx, 1)^2 + ss_locs(choose_idx, 2)^2);
chosen_t = quakes{cell2mat(quakes(:, 1)) == choose_idx, 2}; % day, can be found in "quakes"

figure;
plot(quaket_wk./86400, quaker_wk./1e3, 'r.', 'MarkerSize', 20); hold on;
plot(quaket_as./86400, quaker_as./1e3, 'b+', 'MarkerSize', 10);
plot(chosen_t, chosen_r./1e3, 'r+', 'MarkerSize', 20);
xlabel('days since injection started'); ylabel('distance from injector (km)');
xline(13, 'k--'); yline(0.7, 'k--');


%% plot aseismic and seismic moment after end of injection
addpath(fullfile(above_dir ,'/output/final/'))
load('as_slip_long.mat')
load('M1.mat');
load('M2.mat');
load('Cooper_Basin_Catalog_HAB_4.mat')
event_t = Catalog(2).val; event_t = (event_t-event_t(1)); event_t = event_t-3;
event_M0 = Catalog(7).val;

t_end = M2.H4_time(end)./86400; % time at which injection ends (days)

% compute observed cumulative seismic moment after end of injection----------------
[~, end_idx_ss] = min(abs(event_t - t_end));
postT_ss = event_t(end_idx_ss:end).*86400;
cM_ss = cumsum(event_M0); % measured cumulative seismic moment
post_cM_ss = cM_ss(end_idx_ss:end)-cM_ss(end_idx_ss);

% compute simulated cumulative aseismic moment after end of injection--------------------
mu = M1.mu;
A = M1.dx*M1.dz; % area of one grid element
[~, end_idx_as] = min(abs(t_as_s - t_end*86400));
postT_as = t_as_s(end_idx_as:end);

% displacement field after injection stopped
postD_as = Dx_as_s(:,:,end_idx_as:end); 
postD_as= postD_as - postD_as(:,:,1);
post_cM_as = squeeze(sum(sum(postD_as .* mu .* A, 1), 2));

figure;
plot((postT_as-postT_as(1))./86400, post_cM_as, 'k-', 'LineWidth', 2); hold on;
plot((postT_ss-postT_ss(1))./86400, post_cM_ss, 'k--', 'LineWidth', 2);
xlabel('days after injection stopped'); ylabel('cumulative moment change')
xlim([0, 3]);

