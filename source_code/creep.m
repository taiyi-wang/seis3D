function [t_as_s, psi_as_s, Dx_as_s, p_as_s, taux_as_s, tauy_as_s] = creep(M, time, pressure, pressureflag, plotflag, saveflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run aseismic creeping on main fault using adaptive time stepping

% Input:
% M             = parameter structure
% time          = Nt x 1 time vector associated with pressure diffusion
% pressure      = Nz x Nx x Nt pressure diffusion history on a Cartesian grid
% pressureflag  = 1 for using time-dependent pressure field due to injection, 0 for using constant initial pore pressure
% plotflag     = 1 for plotting slip velocity, shear stress change, pore pressure while computing
% saveflag     = 1 for saving results while running

% Output:
% t_as_s   = a Ntime X 1 vector of time steps (s) from adaptive time stepping
% psi_as_s = Nz x Nx x Ntime array of dimensionless state parameters at the aseismic main fault
% Dx_as_s  = Nz x Nx x Ntime array of cumulative slip at aseismic main fault 
% p_as_s   = Nz x Nx x Ntime array of pore pressure at aseismic main fault
% taux_as_s  = Nz x Nx x Ntime array of shear stress (along dip) at aseismic main fault
% tauy_as_s  = Nz x Nx x Ntime array of normal stress at aseismic main fault
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
disp('computing aseismic slip on main fault');

saveN = 2; % save every 2 steps

M1 = M{1};
M2 = M{2};

tmax = M1.duration;
dt_max = M1.dt_max;
dt_min = M1.dt_min;

tol = M1.tol; % tolerance 
dt = M1.dt;   % initial time step (s)
safety = M1.safetyFactor; % safety factor
     
% Update current values of time, displacement, state, normal stress, pressure,
t_a = 0; D_as_a = M1.Dx0; psi_as_a = M1.s0; 
tauy_as_a = M1.tauy_as_0; 

if pressureflag == 1
    p_a = interp_t(time, pressure, t_a); 
elseif pressureflag == 0
    p_a = reshape(M2.p_pore, M2.mz, M2.mx); 
end

% Stage 1 update
[V1_as, G1_as, tau1_x_as, tau1_y_as] = stateODE(D_as_a, psi_as_a, t_a, p_a, tauy_as_a, M1);
    
% Store variables
t_as_s = t_a;
Vx_as_s(:,:,1) = V1_as(:,:);
taux_as_s(:,:,1) = tau1_x_as(:,:);
tauy_as_s(:,:,1) = tau1_y_as(:,:);
Dx_as_s(:,:,1) = D_as_a(:,:);
psi_as_s(:,:,1) = psi_as_a;
p_as_s(:,:, 1) = p_a;

% current error and time step
err_a=0; dt_a=dt;
err_s=err_a; dt_s=dt_a;

while t_a < tmax
    % adjust dt to stop at tmax
    if t_a+dt_a>tmax, dt_a=tmax-t_a; end
        
    % three-stage method with embedded error estimate
    % Note: 
    % p_a needs to be updated at intermediate times (currently not implemented, because inaccuracy in time interpolation causes numerical error accumulation)
    % tauy_as_a currently is fixed to initial values, so does not need to be updated at intermediate times
    %if pressureflag == 1
        %p_a = interp_t(time, pressure, t_a+0.5*dt_a);
    %end
    [V2_as, G2_as] = stateODE(D_as_a+0.5*dt_a*V1_as, psi_as_a+0.5*dt_a*G1_as, t_a+0.5*dt_a, p_a, tauy_as_a, M1);
    
    %if pressureflag == 1
        %p_a = interp_t(time, pressure, t_a+dt_a); 
    %end
    [V3_as, G3_as] = stateODE(D_as_a+dt_a*(-V1_as+2*V2_as), psi_as_a+dt_a*(-G1_as+2*G2_as), t_a+dt_a, p_a, tauy_as_a, M1);
        
    % second order update
    D2_as   = D_as_a+dt_a/2*(V1_as+V3_as);
    psi2_as = psi_as_a+dt_a/2*(G1_as+G3_as);
        
    % third order update
    D3_as   = D_as_a+dt_a/6*(V1_as+4*V2_as+V3_as);
    psi3_as = psi_as_a+dt_a/6*(G1_as+4*G2_as+G3_as);
        
    q = 2; % order of accuracy of lower order update
        
    % local error estimate
    err_a = norm([D2_as(:)-D3_as(:); psi2_as(:)-psi3_as(:)])./sqrt(M1.Nx.*M1.Nz);
        
    if err_a<tol
        % use third-order update for the solution
        t_a = t_a+dt_a; D_as_a = D3_as; psi_as_a = psi3_as; 
        
        % store solution
        t_as_s=[t_as_s; t_a];
        Dx_as_s(:,:,length(t_as_s)) = D_as_a(:,:);
        psi_as_s(:,:,length(t_as_s)) = psi_as_a(:,:); 
        
        % store error and time step
        err_s=[err_s; err_a]; dt_s=[dt_s; dt_a];
            
        % evaluate stage 1 values for next time step
        [V1_as, G1_as, tau1_x_as, tau1_y_as] = stateODE(D_as_a, psi_as_a, t_a, p_a, tauy_as_a, M1);
        
        % stage 1 values are stored
        Vx_as_s(:,:,length(t_as_s)) = V1_as(:,:); 
        taux_as_s(:,:,length(t_as_s)) = tau1_x_as(:,:); tauy_as_s(:,:,length(t_as_s)) = tau1_y_as(:,:);
        
        % get pressure for current time
        if pressureflag == 1
            p_a = interp_t(time, pressure, t_a); 
        end
        p_as_s(:,:,length(t_as_s)) = p_a;
            
        % save results to a mat file periodically
        stepN = length(t_as_s);
        if (saveflag == 1) && ((stepN/saveN) == round(stepN/saveN))
            if isfile('as_output.mat')
                output = matfile('as_output.mat','Writable',true);          
                output.t_as_s((stepN-saveN+1):stepN, 1) = t_as_s((stepN-saveN+1):stepN, 1);  
                output.err_s((stepN-saveN+1):stepN, 1) = err_s((stepN-saveN+1):stepN, 1);
                output.dt_s((stepN-saveN+1):stepN, 1) = dt_s((stepN-saveN+1):stepN, 1);
                output.psi_as_s(:,:,(stepN-saveN+1):stepN) = psi_as_s(:,:,(stepN-saveN+1):stepN); 
                output.Dx_as_s(:,:,(stepN-saveN+1):stepN) = Dx_as_s(:,:,(stepN-saveN+1):stepN);
                output.Vx_as_s(:,:,(stepN-saveN+1):stepN) = Vx_as_s(:,:,(stepN-saveN+1):stepN);
                output.taux_as_s(:,:,(stepN-saveN+1):stepN) = taux_as_s(:,:,(stepN-saveN+1):stepN);
                output.tauy_as_s(:,:,(stepN-saveN+1):stepN) = tauy_as_s(:,:,(stepN-saveN+1):stepN);
                output.p_as_s(:,:,(stepN-saveN+1):stepN) = p_as_s(:,:,(stepN-saveN+1):stepN);
            else
                save('as_output.mat', 't_as_s', 'psi_as_s', 'Dx_as_s', 'Vx_as_s', 'taux_as_s', 'tauy_as_s', 'err_s', 'dt_s', 'p_as_s', '-v7.3'); % create an object - output
            end
        end
        
        if plotflag == 1
            h1 = subplot(1,3,1);
            surf(M1.X./1000, M1.Z./1000, log10(V1_as),'EdgeColor', 'none');%caxis('manual');caxis([10e-10, 10e1]);
            xlim([min(M1.X(:)./1000), max(M1.X(:)./1000)]);ylim([min(-M1.Z(:)./1000), max(-M1.Z(:)./1000)]);
            colorbar;view([0, 90]);colormap(h1, 'hot');
            title(sprintf('$log(V_x)$ at t = %f days', t_a/86400), 'Interpreter', 'latex');
            xlabel('along dip (km)');ylabel('along strike (km)')
            
            h2 = subplot(1,3,2);
            surf(M1.X./1000, M1.Z./1000, (tau1_x_as - squeeze(taux_as_s(:,:,1)))./1e6,'EdgeColor', 'none');%caxis('manual');caxis([10e-10, 10e1]);
            xlim([min(M1.X(:)./1000), max(M1.X(:)./1000)]);ylim([min(-M1.Z(:)./1000), max(-M1.Z(:)./1000)]);
            colorbar;view([0, 90]);colormap(h2, 'winter');
            title(sprintf('$\\Delta \\tau_{yx}$ (MPa) at t = %f days', t_a/86400), 'Interpreter', 'latex');
            xlabel('along dip (km)');ylabel('along strike (km)')
            
            h3 = subplot(1,3,3);
            surf(M1.X./1000, M1.Z./1000,(p_a-p_as_s(:,:,1))./1e6,'EdgeColor', 'none');
            xlim([min(M1.X(:)./1000), max(M1.X(:)./1000)]);ylim([min(-M1.Z(:)./1000), max(-M1.Z(:)./1000)]);
            colorbar;view([0, 90]);colormap(h3, 'summer');
            title(sprintf('$\\Delta p$ (MPa) at t = %f days', t_a/86400), 'Interpreter', 'latex');
            xlabel('along dip (km)');ylabel('along strike (km)')
            
            drawnow
        end    
    end
    
    % adjust time step
    dt_a = safety*dt_a*(tol/err_a)^(1/(q+1));
    dt_a = min(dt_a,dt_max);
    
end

disp('finished')
end
    




