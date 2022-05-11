% Make plots for publication
% Taiyi Wang 11/30/2021

current_dir  = pwd;
idcs   = strfind(current_dir ,'/');
above_dir = current_dir(1:idcs(end)-1); % one level above current directory

%addpath(fullfile(above_dir, '/output'))
addpath(fullfile(above_dir ,'/source_code'))
addpath(fullfile(above_dir ,'/input'))

%% Adjust font and line
set(0,'DefaultTextFontSize',16)
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultLineLineWidth', 2);
set(0,'DefaultLineMarkerSize',12);

%% Load data needed
load('M1.mat'); 
load('M2.mat');
load('injection_output_15')
load('ss_output');
load('as_output');
%% Make schematic to show the main fault and seismic patches
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


figure;
for i = 1:2:length(new_t)
    plot(Z(:, xx_0_idx)'./1000, new_dp_x(:,i)./1e6, 'r', 'LineWidth', 1);hold on;  
    plot(X(zz_0_idx,:)./1000, new_dp_z(:,i)./1e6, 'r--', 'LineWidth', 1)
end
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

%% Plot seismic slip along each dimension
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
% Unpack variables
Nx = M1.Nx; Nz = M1.Nz; p_pore = M2.p_pore;
N_ss = M1.N_ss; mu = M1.mu;

mid_x_idx = (Nx+1)/2;
mid_z_idx = (Nz+1)/2;

% Velocity on main fault in log10(m/s)
Vx_as = squeeze(Vx_as_s(:,mid_z_idx,:));
Vx_as_max = max(Vx_as, [],'all');
Vx_as_min = min(Vx_as, [],'all');

z_locs = Z(:, 1);
xlocs_ss = ss_locs(:,1)./1000; % x coordinates of the spring sliders in km
zlocs_ss = ss_locs(:,2)./1000; % x coordinates of the spring sliders in km
[T_grid, Z_grid] = meshgrid(t_as_s, z_locs);

% retrieve locations of observed seismic events along z = 0 km, for x >= 0 km
load('Cooper_Basin_Catalog_HAB_4.mat')
event_lat = Catalog(4).val; event_lon = Catalog(5).val; elevation = Catalog(6).val; event_t = Catalog(2).val;
event_t = (event_t-event_t(1)); event_Mw = Catalog(8).val; event_M0 = Catalog(7).val;
HB4_lat = -27.8115; HB4_lon = 140.7596;% from well completion report

% convert event locations to meters with regard to HB4
event_local = llh2local([event_lon'; event_lat'; elevation'],[HB4_lon; HB4_lat; 0]);
event_x = event_local(1,:).*1000; event_z = event_local(2,:).*1000; 

% find events along x-axis (west-east), extending from HB4 towards west-east
event_narroaw_idx = intersect(find(event_z<100 & event_z>-100), find(event_Mw > 0));

figure;
p = pcolor(T_grid./86400, Z_grid./1000, log10(Vx_as)); colormap(parula);colorbar; caxis([log10(Vx_as_min), log10(Vx_as_max)-2]);
p.FaceColor = 'interp'; set(p, 'EdgeColor', 'none'); p.FaceAlpha = 1; hold on;

for i = 1:length(event_narroaw_idx)
    plt = plot(event_t(event_narroaw_idx(i))-3, event_x(event_narroaw_idx(i))./1000, 'w.', 'MarkerSize', event_Mw(event_narroaw_idx(i))*10);% seismic catalogue starts on 11/10/2012, whereas pressure data starts on 11/13/2012
    plt.Color(4) = 0.2; 
end

quakes = find_quakes(t_ss_s, Vx_ss_s, Dx_ss_s, M1);

for i = 1:size(quakes, 1)
    plot(quakes{i, 2}, xlocs_ss(i), '.', 'Color', [0.8500, 0.3250, 0.0980, 0.2], 'MarkerSize', abs(quakes{i, 4}).*10);  
end

% make a legend for seismic event moment magnitude
plot(1.1, -0.6, 'w.', 'MarkerSize', 2*10); %text(1.2, -0.5, 'M_w = 2', 'Color', 'w');
plot(1.1, -0.75, 'w.', 'MarkerSize', 1.5*10); %text(1.2, -0.65, 'M_w = 1.5', 'Color', 'w');
plot(1.1, -0.9, 'w.', 'MarkerSize', 1*10); %text(1.2, -0.8, 'M_w = 1', 'Color', 'w');
view([0, 90]);
xlim([min(t_s./86400), 17]); ylim([-1, 1]);
ylabel('z (km)'); xlabel('time (days)'); %title('Color scale = aseismic slip velocity (cm/day); contour = perturbation pressure (MPa)')

%% Histograms of seismicity
% 1. time dependence of seismicity
data_event_t = event_t-3; % seismic catalogue starts on 11/10/2012, whereas pressure data starts on 11/13/2012
data_event_t = data_event_t(data_event_t<=17);
ss_event_t = cell2mat(quakes(:,2));

figure;
histogram(data_event_t, 'Normalization', 'Probability', 'FaceAlpha', 1, 'FaceColor', [1, 1, 1], 'LineWidth', 1, 'BinWidth', 0.2); hold on;
histogram(ss_event_t, 'Normalization', 'Probability', 'LineWidth', 1, 'FaceAlpha', 0.8, 'BinWidth', 0.2);
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
histogram(event_x(data_narroaw_z_idx)./1000, 'Normalization', 'Probability', 'FaceAlpha', 1, 'FaceColor', [1, 1, 1], 'LineWidth', 1, 'BinWidth', 0.03);hold on;
histogram(pred_x(pred_narroaw_z_idx)./1000, 30, 'Normalization', 'Probability', 'LineWidth', 1, 'FaceAlpha', 0.8, 'BinWidth', 0.03); 
xlim([-1.5, 1.5]);
xlabel('distance along dip (km)');
ylabel('normalized frequency')

figure;
histogram(event_z(data_narroaw_x_idx)./1000, 'Normalization', 'Probability', 'FaceAlpha', 1, 'FaceColor', [1, 1, 1], 'LineWidth', 1, 'BinWidth', 0.03);hold on;
histogram(pred_z(pred_narroaw_x_idx)./1000, 30, 'Normalization', 'Probability', 'LineWidth', 1, 'FaceAlpha', 0.8, 'BinWidth', 0.03); 
xlim([-1.5, 1.5]);
xlabel('distance along strike (km)');
ylabel('normalized frequency')

%% Normalized cumulative seismic/aseismic moment release
A = M1.dx .* M1.dz;
mu = M1.mu;

% 1. compare data/prediction in seismic moment release
t_m = (0:0.05:17).*86400; % evenly spaced time for resampling

% compute predicted total seismic moment release (default case)
cM_ss_p = zeros(1, length(t_m));
load('ss_output_15.mat', 't_ss_s', 'Dx_ss_s')
for i = 1:N_ss
    t = t_ss_s{i};
    d = Dx_ss_s{i};
    m = A .* mu .*  d; % cumulative moment
    
    cM_ss_p = cM_ss_p + interp1(t, m, t_m);
end
clear t d m t_ss_s Dx_ss_s

% compute predicted total seismic moment release (12 MPa initial shear stress)
cM_ss_p_12 = zeros(1, length(t_m));
load('ss_output_12.mat', 't_ss_s', 'Dx_ss_s')
for i = 1:N_ss
    t = t_ss_s{i};
    d = Dx_ss_s{i};
    m = A .* mu .*  d; % cumulative moment
    
    cM_ss_p_12 = cM_ss_p_12 + interp1(t, m, t_m);
end
clear t d m t_ss_s Dx_ss_s

% compute predicted total seismic moment release (17 MPa initial shear stress)
cM_ss_p_17 = zeros(1, length(t_m));
load('ss_output_17.mat', 't_ss_s', 'Dx_ss_s')
for i = 1:N_ss
    t = t_ss_s{i};
    d = Dx_ss_s{i};
    m = A .* mu .*  d; % cumulative moment
    
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


%% Plot the 2D distribution of seismicity at any given time
chosen_day = 10;

figure;

event_idx = find(event_t-3 -chosen_day <= 1 & event_t-3 -chosen_day >= 0);
plot(event_x(event_idx)./1000, event_z(event_idx)./1000, 'k.');
[~, taux_day_idx] = min(abs((t_as_s./86400)-chosen_day));
dtaux_day = squeeze(taux_as_s(:,:,taux_day_idx)) - squeeze(taux_as_s(:,:,1));
h = pcolor(X./1000, Z./1000, dtaux_day./1e6);
set(h, 'EdgeColor', 'none');
hold on;

for i = 1:N_ss
    Nevent_day = find(pop_list{i, 1}-chosen_day <= 1 & pop_list{i, 1}-chosen_day >= 0); % get number of events in the chosen day
    if Nevent_day >= 1
        plot(zlocs_ss(i), xlocs_ss(i), 'r.');
    end
end
xlabel('along dip (SSW-NNE)');ylabel('along strike')
colorbar; 
xlim([-2, 2]); ylim([-2, 2])
clear i event_idx Nevent_day

%% spring slider phase diagram
% choose the index of spring slider
choose_idc = 87;

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
    
    semilogy(t_ss_s{choose_idx}, Vx_ss_s{choose_idx}, 'k-'); hold on;
end

%% spring slider pressure and aseismic loading stress history
choose_idx = 87;

% compute total strength change
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
plot(t_as_s./86400, M1.s0_ss.*normal_stress_change(2:end)./1e6, 'b'); % shear strength change due to normal stress change produced by aseismic loading
plot(t_ss_s{choose_idx}./86400, -M1.s0_ss.*(p_ss_s{choose_idx} - p_ss_s{choose_idx}(1))./1e6, 'r'); % shear strength change due to normal stress change produced by pressure diffusion
legend('\delta \tau_{str}', '\delta \tau_{as}', '\Psi_0 \cdot \delta \sigma_{as}', '-\Psi_0 \cdot \delta p_{ss}')
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

%%
figure;
scatter3(event_x./1e3, event_z./1e3, elevation', (event_Mw./max(event_Mw)+0.1).*5e2, event_t, '.');hold on; colorbar; 
plot3(repelem(0, 5, 1), repelem(0, 5, 1), linspace(-3.5, -4.085, 5), 'Color', [0.2, 0.4, 0.8], 'LineWidth', 2);
xlabel('East (km)'); ylabel('North (km)');
colormap('jet'); axis equal;
grid on;



