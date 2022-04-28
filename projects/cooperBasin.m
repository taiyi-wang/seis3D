%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script sets up a simulation for the spatial-temporal distributions of 
% Cooper Basin Enhanced Geothermal System (EGS) induced seismicity for the
% injection period in 2012 at Habanero 4 well

% Taiyi Wang 06/29/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add paths
current_dir  = pwd;
idcs   = strfind(current_dir ,'/');
above_dir = current_dir(1:idcs(end)-1); % one level above current directory

addpath(fullfile(above_dir, 'source_code/.'))
addpath(fullfile(above_dir, 'input/.'));
addpath(fullfile(above_dir, 'output/.'));
%% Options
loadflag = 0; % load a previous setup
testflag = 0; % not benchmark testing
rglgridflag = 0; % use logrithmic radial grid
plotflag = 1; % plotting in real time 
saveflag = 1; % save the main fault slip simulation results

%% Set up parameters
if loadflag == 1
    load('M1.mat'); load('M2.mat')
else
    [M1, M2] = setup_model(); % model parameters
    M2 = setup_data(M2, above_dir);
end

%% Save the input parameters
if ~loadflag  % save input if setting up in real time
    save('M1.mat', 'M1')
    save('M2.mat', 'M2')
end
%% Run simulation
[r, t_s, p_s, pWH, dpw0, q_s] = injection2(M2, testflag, rglgridflag, plotflag, saveflag); % linear pressure diffusion simulation 
%%
p_s_xy = lograd_to_rgxy(r, t_s, p_s, M1, saveflag);                                                  % interpolate pressure onto regular, Cartesian grid
%%
[t_as_s, psi_as_s, Dx_as_s, p_as_s, taux_as_s, tauy_as_s] = creep({M1, M2}, t_s', p_s_xy, 1, plotflag, saveflag);    % aseismic creep simulation
%%
[t_ss_s, psi_ss_s, Vx_ss_s, Dx_ss_s, p_ss_s, taux_ss_s, idx_timeout] = seismicity(M1, M2, t_as_s, Dx_as_s, p_as_s, saveflag); % off-fault seismicity simulation



